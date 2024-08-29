__version__ = "0.2.26"
# Data version is kept in generate_transcript_version.json_schema_version

def get_data_schema_int(version: str) -> int:
    """ Return an int which increments upon breaking changes - ie anything other than patch """
    major, minor, patch = version.split(".")
    return 1000 * int(major) + int(minor)
