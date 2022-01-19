__version__ = "0.1.1"


def get_json_schema_version():
    """ Return an int which increments upon breaking changes - ie anything other than patch """
    major, minor, patch = __version__.split(".")
    return 1000 * int(major) + int(minor)
