import unittest
import itertools

from zigzag.sequences import *

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1e6)

class TestSequences(unittest.TestCase):
	"""Test sequence manipulation functions."""
	
	def _generate_context(self, width):
		"""Generate context for sequence alignment."""
		
		half_width = width // 2
		segment_length_cases = {
			0,
			half_width - 1,
			half_width,
			half_width + 1,
		}
		pattern_length_cases = {
			1,
			segment_length - 1,
			segment_length,
		}
		ambiguity_cases = (True, False)


		return context

	def _generate_patterns(self):


		possible_characters = ['N', 'C', 'P', 'X', '-']
		

		pass

	def _generate_expected_sequences(self):
		pass

	def test_align_sequences_even(self):
		"""Test align_sequences given even aligned sequence width."""

		# Arbitrary width. Using 8 based on MEROPS substrate width.
		width = 8
		

	def test_align_sequences_odd(self):
		"""Test align_sequences given odd aligned sequence width."""

		# Arbitrary width. Using 8 based on Motif-X default width.
		width = 13
