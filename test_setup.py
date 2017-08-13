from mock import patch, MagicMock

from scripts import base_setup


def test_get_package_data():
    print(base_setup.get_package_data())

def test_get_version_info():
    print(base_setup.get_version_info())


def test_get_pyx_pxd():
    print(base_setup.get_pyx_pxd())
