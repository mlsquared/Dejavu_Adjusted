import multiprocessing
import os
import sys
import traceback
from itertools import groupby
from time import time
from typing import Dict, List, Tuple

import dejavu.logic.decoder as decoder
from dejavu.base_classes.base_database import get_database
from dejavu.config.settings import (DEFAULT_FS, DEFAULT_OVERLAP_RATIO,
                                    DEFAULT_WINDOW_SIZE, FIELD_FILE_SHA1,
                                    FIELD_TOTAL_HASHES,
                                    FINGERPRINTED_CONFIDENCE,
                                    FINGERPRINTED_HASHES, HASHES_MATCHED,
                                    INPUT_CONFIDENCE, INPUT_HASHES, OFFSET,
                                    OFFSET_SECS, SONG_ID, SONG_NAME, TOPN)
from dejavu.logic.fingerprint import fingerprint

from .fingerprint_writer_custom import FingerprintWriter


class Dejavu:
    def __init__(self, config):
        self.config = config

        # initialize db
        db_cls = get_database(config.get("database_type", "mysql").lower())

        self.db = db_cls(**config.get("database", {}))
        self.db.setup()

        # if we should limit seconds fingerprinted,
        # None|-1 means use entire track
        self.limit = self.config.get("fingerprint_limit", None)
        if self.limit == -1:  # for JSON compatibility
            self.limit = None
        self.__load_fingerprinted_audio_hashes()

    def __load_fingerprinted_audio_hashes(self) -> None:
        """
        Keeps a dictionary with the hashes of the fingerprinted songs, in that way is possible to check
        whether or not an audio file was already processed.
        """
        # get songs previously indexed
        self.songs = self.db.get_songs()
        self.songhashes_set = set()  # to know which ones we've computed before
        for song in self.songs:
            song_hash = song[FIELD_FILE_SHA1]
            self.songhashes_set.add(song_hash)

    def get_fingerprinted_songs(self) -> List[Dict[str, any]]:
        """
        To pull all fingerprinted songs from the database.

        :return: a list of fingerprinted audios from the database.
        """
        return self.db.get_songs()

    def delete_songs_by_id(self, song_ids: List[int]) -> None:
        """
        Deletes all audios given their ids.

        :param song_ids: song ids to delete from the database.
        """
        self.db.delete_songs_by_id(song_ids)

    def fingerprint_directory(self, path: str, extensions: str, nprocesses: int = None) -> None:
        """
        Given a directory and a set of extensions it fingerprints all files that match each extension specified.

        :param path: path to the directory.
        :param extensions: list of file extensions to consider.
        :param nprocesses: amount of processes to fingerprint the files within the directory.
        """
        # Try to use the maximum amount of processes if not given.
        try:
            nprocesses = nprocesses or multiprocessing.cpu_count()
        except NotImplementedError:
            nprocesses = 1
        else:
            nprocesses = 1 if nprocesses <= 0 else nprocesses

        pool = multiprocessing.Pool(nprocesses)

        filenames_to_fingerprint = []
        for filename, _ in decoder.find_files(path, extensions):
            # don't refingerprint already fingerprinted files
            if decoder.unique_hash(filename) in self.songhashes_set:
                print(f"{filename} already fingerprinted, continuing...")
                continue

            filenames_to_fingerprint.append(filename)

        # Prepare _fingerprint_worker input
        worker_input = list(zip(filenames_to_fingerprint, [self.limit] * len(filenames_to_fingerprint)))

        # Send off our tasks
        iterator = pool.imap_unordered(Dejavu._fingerprint_worker, worker_input)

        # Loop till we have all of them
        while True:
            try:
                song_name, hashes, file_hash, peaks = next(iterator)
            except multiprocessing.TimeoutError:
                continue
            except StopIteration:
                break
            except Exception:
                # Print traceback because we can't reraise it here
                traceback.print_exc(file=sys.stdout)
            else:
                # Write to .fpt file instead of adding to DB
                self.custom_write_fpt_file(song_name, hashes, file_hash, peaks)

                # sid = self.db.insert_song(song_name, file_hash, len(hashes))
                #
                # self.db.insert_hashes(sid, hashes)
                # self.db.set_song_fingerprinted(sid)
                # self.__load_fingerprinted_audio_hashes()

        pool.close()
        pool.join()

    @staticmethod
    def custom_write_fpt_file(filename, hashes, file_hash, peaks):
        print("Now in custom function")
        working_dir = "/" + os.getcwd().split("/")[1]  # Should be "/app"
        fpt_dir = os.getenv("FPT_FILES_DIRECTORY")
        fpt_dir = os.path.join(working_dir, fpt_dir)

        # Sort the hashes based on offset first, then hash value
        hashes = sorted(hashes, key=lambda x: (x[1], x[0]))


        fields = [("hash", "str"), ("ts", "int"), ("features", "str")]
        last_offset = 0 if len(hashes) == 0 else hashes[-1][1]
        metadata = {
            "algorithm": "dejavu",
            "content_id": filename,
            "audio_duration": round(float(last_offset) / DEFAULT_FS * DEFAULT_WINDOW_SIZE * DEFAULT_OVERLAP_RATIO, 3),
            "num_fingerprints": len(hashes),
            "filehash": file_hash
        }

        writer = FingerprintWriter(filename=os.path.join(fpt_dir, filename + ".fpt"), fields=fields, kwargs=metadata)
        print(f"Created fpt file, {os.path.join(fpt_dir, filename + '.fpt')}")
        # Adding some comments first to the .fpt file
        # These are unique to dejavu algorithm, change as needed
        writer.add_comments("")
        writer.add_comments("Duration is a float measured in seconds of the **APPROXIMATED** audio duration")
        writer.add_comments("Unless the ending is Silent/Near Silent (current guess), Duration should be very close approximation of actual duration")
        writer.add_comments("")
        writer.add_comments("Timestamp (ts) is measured in milliseconds")
        writer.add_comments("")
        writer.add_comments("Features are: freq1, freq2, t_delta")
        writer.add_comments("")


        for hash_val, offset in hashes:
            # Convert the offset into milliseconds
            ts = round(float(offset) / DEFAULT_FS * DEFAULT_WINDOW_SIZE * DEFAULT_OVERLAP_RATIO * 1000)
            hash_key = (hash_val, offset)

            # If Hash isn't present in peaks, write (NULL) to features, this should not happen.
            if hash_key not in peaks:
                hash_val = "0x" + hash_val.upper()  # Adjust the hash
                writer.write_fingerprints(hash_val, ts, "(NULL)")
                continue

            hash_val = "0x" + hash_val.upper()  # Adjust the hash
            peak_entries = peaks[hash_key]

            # Convert the peak values to a string representation
            peak_values = "(" + ",".join(str(int(a)) for a in peak_entries[0]) + ")"

            writer.write_fingerprints(hash_val, ts, peak_values)

    def fingerprint_file(self, file_path: str, song_name: str = None) -> None:
        """
        Given a path to a file the method generates hashes for it and stores them in the database
        for later be queried.

        :param file_path: path to the file.
        :param song_name: song name associated to the audio file.
        """
        song_name_from_path = decoder.get_audio_name_from_path(file_path)
        song_hash = decoder.unique_hash(file_path)
        song_name = song_name or song_name_from_path
        # don't refingerprint already fingerprinted files
        if song_hash in self.songhashes_set:
            print(f"{song_name} already fingerprinted, continuing...")
        else:
            song_name, hashes, file_hash, peaks = Dejavu._fingerprint_worker((file_path, self.limit)) # Modify pass tuple
            self.custom_write_fpt_file(song_name, hashes, file_hash, peaks)

            # sid = self.db.insert_song(song_name, file_hash)

            # self.db.insert_hashes(sid, hashes)
            # self.db.set_song_fingerprinted(sid)
            # self.__load_fingerprinted_audio_hashes()

    def generate_fingerprints(self, samples: List[int], Fs=DEFAULT_FS) -> Tuple[List[Tuple[str, int]], float]:
        f"""
        Generate the fingerprints for the given sample data (channel).

        :param samples: list of ints which represents the channel info of the given audio file.
        :param Fs: sampling rate which defaults to {DEFAULT_FS}.
        :return: a list of tuples for hash and its corresponding offset, together with the generation time.
        """
        t = time()
        hashes = fingerprint(samples, Fs=Fs)
        fingerprint_time = time() - t
        return hashes, fingerprint_time

    def find_matches(self, hashes: List[Tuple[str, int]]) -> Tuple[List[Tuple[int, int]], Dict[str, int], float]:
        """
        Finds the corresponding matches on the fingerprinted audios for the given hashes.

        :param hashes: list of tuples for hashes and their corresponding offsets
        :return: a tuple containing the matches found against the db, a dictionary which counts the different
         hashes matched for each song (with the song id as key), and the time that the query took.

        """
        t = time()
        matches, dedup_hashes = self.db.return_matches(hashes)
        query_time = time() - t

        return matches, dedup_hashes, query_time

    def align_matches(self, matches: List[Tuple[int, int]], dedup_hashes: Dict[str, int], queried_hashes: int,
                      topn: int = TOPN) -> List[Dict[str, any]]:
        """
        Finds hash matches that align in time with other matches and finds
        consensus about which hashes are "true" signal from the audio.

        :param matches: matches from the database
        :param dedup_hashes: dictionary containing the hashes matched without duplicates for each song
        (key is the song id).
        :param queried_hashes: amount of hashes sent for matching against the db
        :param topn: number of results being returned back.
        :return: a list of dictionaries (based on topn) with match information.
        """
        # count offset occurrences per song and keep only the maximum ones.
        sorted_matches = sorted(matches, key=lambda m: (m[0], m[1]))
        counts = [(*key, len(list(group))) for key, group in groupby(sorted_matches, key=lambda m: (m[0], m[1]))]
        songs_matches = sorted(
            [max(list(group), key=lambda g: g[2]) for key, group in groupby(counts, key=lambda count: count[0])],
            key=lambda count: count[2], reverse=True
        )

        songs_result = []
        for song_id, offset, _ in songs_matches[0:topn]:  # consider topn elements in the result
            song = self.db.get_song_by_id(song_id)

            song_name = song.get(SONG_NAME, None)
            song_hashes = song.get(FIELD_TOTAL_HASHES, None)
            nseconds = round(float(offset) / DEFAULT_FS * DEFAULT_WINDOW_SIZE * DEFAULT_OVERLAP_RATIO, 5)
            hashes_matched = dedup_hashes[song_id]

            song = {
                SONG_ID: song_id,
                SONG_NAME: song_name.encode("utf8"),
                INPUT_HASHES: queried_hashes,
                FINGERPRINTED_HASHES: song_hashes,
                HASHES_MATCHED: hashes_matched,
                # Percentage regarding hashes matched vs hashes from the input.
                INPUT_CONFIDENCE: round(hashes_matched / queried_hashes, 2),
                # Percentage regarding hashes matched vs hashes fingerprinted in the db.
                FINGERPRINTED_CONFIDENCE: round(hashes_matched / song_hashes, 2),
                OFFSET: offset,
                OFFSET_SECS: nseconds,
                FIELD_FILE_SHA1: song.get(FIELD_FILE_SHA1, None).encode("utf8")
            }

            songs_result.append(song)

        return songs_result

    def recognize(self, recognizer, *options, **kwoptions) -> Dict[str, any]:
        r = recognizer(self)
        return r.recognize(*options, **kwoptions)

    @staticmethod
    def _fingerprint_worker(arguments):
        # Pool.imap sends arguments as tuples so we have to unpack
        # them ourself.
        try:
            file_name, limit = arguments
        except ValueError:
            pass

        song_name, extension = os.path.splitext(os.path.basename(file_name))

        fingerprints, peaks, file_hash = Dejavu.get_file_fingerprints(file_name, limit, print_output=True)

        return song_name, fingerprints, file_hash, peaks

    @staticmethod
    def get_file_fingerprints(file_name: str, limit: int, print_output: bool = False):
        channels, fs, file_hash = decoder.read(file_name, limit)
        fingerprints = set()
        peaks = {}
        channel_amount = len(channels)
        for channeln, channel in enumerate(channels, start=1):
            if print_output:
                print(f"Fingerprinting channel {channeln}/{channel_amount} for {file_name}")

            peaks_hashes = fingerprint(channel, Fs=fs)

            if print_output:
                print(f"Finished channel {channeln}/{channel_amount} for {file_name}")

            channel_peaks, hashes = peaks_hashes[1], peaks_hashes[0]

            fingerprints |= set(hashes)
            
            for key, value in channel_peaks.items():
                if key not in peaks:
                    peaks[key] = value
                else:
                    peaks[key].extend(value)

        return fingerprints, peaks, file_hash
