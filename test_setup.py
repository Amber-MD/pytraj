from mock import patch, MagicMock

from scripts import base_setup


def test_get_package_data():
    print(base_setup.get_package_data())
