from django.test import TestCase

class FileUploadTestCase(TestCase):
	"""Test file functionality."""

	def _create_test_file(self, path):
		"""Generate a one-line text file for testing."""
		with open(path, 'w') as f:
			f.write('ABDERKVI\n')

	def _delete_test_file(self, path):
		"""Delete test file."""
		pass

	def test_file_upload(self):
		test_file_path  = '/tmp/test_file.txt'
		data = self._create_test_file(test_file_path)