import json
import os


class FingerprintWriter:
    """
    This is my implementation based on the `fingerprint_writer.py` file in `https://github.com/mlsquared/tools_rd/blob/main/common/`

    This class
    'Writes fingerprints, their timestamps, and associated fields to a fingerprint file.
    The output format is designed to be read by FingerprintReader.'

    Staying faithful to the original file, the default fields would also be

        fingerprints, ts, content_id

    Can be adjusted as needed by passing in a `fields` list. The expected format would be a list of tuples, each
    tuple being (field name, type) pairs.
    """

    def __init__(self, filename: str, fields=None, overwrite_files=False, **kwargs):
        """
        Constructs a FingerprintWriter object that writes to specified file.

        :param filename: A string representing the output file name 
        :param fields: An optional parameter as described above
        :param overwrite_files: A boolean that raises error if given filename already exists. Gives warnings if True
        :param kwargs: Other metadata to be included in the commented JSON dict at the top of the .fpt file
        """

        if fields is None:
            fields = [("fingerprint", "fingerprint"), ("ts", "int"), ("content_id", "str")]
        self.fields = fields

        self.filename = self.get_filename(filename)
        self.wrote_fingerprints = False

        if os.path.exists(self.filename) and not overwrite_files:
            raise FileExistsError(f"{self.filename} already exists. Use a different filename or remove the existing file. Use overwrite_files=True if desired.")

        header_dict = kwargs["kwargs"]
        header_dict["fields"] = self.fields
        header_str = json.dumps(header_dict, indent=4)
        header_str = "# " + ("\n# ".join(header_str.split("\n"))) + "\n"
        header_str += "#\n# (blank comment line above ends the header JSON dict)\n"

        with open(self.filename, "w") as f:
            f.write(header_str)

    @staticmethod
    def get_filename(filename: any) -> str:
        """
        Verifies if the given filename is valid.

        :param filename: Expects a string representing the output file name
        :return: A string for filename that has been adjusted if needed
        """
        if not isinstance(filename, str):
            raise TypeError(f"Invalid filename. Must be a string, {filename=}")
        if filename is None:
            raise ValueError("Must provide a filename")
        if not filename[-1].isalnum() or len(filename) < 1:
            raise ValueError("Filename cannot be empty and last character must be alphanumeric")

        split = filename.split("/")[-1].split(".")
        if len(split) == 1:
            return filename + ".fpt"

        # Not sure if .fpt extension should be enforced. Will comment it out for now.
        # if split[-1] != "fpt":
        #     return "".join(split[:-1]) + ".fpt"

        return filename


    def write_fingerprints(self, *args) -> None:
        """
        Writes the given fingerprint to the associated file given *args

        :param args: An iterable of length equal to fields. Each element in *args should match to the corresponding field
        :return: None
        """

        assert len(args) == len(self.fields), f"Length args={len(args)} != Length fields={len(self.fields)}"

        with open(self.filename, "a") as f:
            try:
                hash_value = str(args[0])
                time = str(args[1])
                features = args[2]
                f.write(f'{hash_value} {time} {features}\n')
                self.wrote_fingerprints = True  # We are now writing fingerprints. Should not be able to add comments now
            except Exception as e:
                raise ValueError(f"Error converting arguments to strings: {e}")


    def add_comments(self, comment: str) -> None:
        """
        Used to add a comment at the head of the .fpt file

        :param comment: The string comment that should be added
        :return: None
        """

        assert isinstance(comment, str), "Given comment should be of string dtype!"
        assert not self.wrote_fingerprints, "Cannot add comments once fingerprints are added.\nCall this method if wish to add comments to header."

        with open(self.filename, "a") as f:
            f.write("# " + comment + "\n")

 