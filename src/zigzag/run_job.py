def run_job(job_id):
	"""Initiate pattern extraction workflow using Django RQ."""
	# Get job row from Jobs table using job_id.
	# Import context from uploaded file or load default Swiss-Prot FASTA.
	# Instantiate Background from context sequences.
	# Import input data from uploaded file.
	# Instantiate PatternContainer with input data and background.
	# Run pattern post-processing.
	# Retrieve (or generate) protease patterns.
	# Generate protease labels.
	# Generate protease pattern frequency matrix.
	# Generate heatmap.
	# Write output file paths to Results table.
	# Report job complete.