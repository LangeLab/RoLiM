import os
import uuid

from django.db import models
from django.dispatch import receiver
from django.contrib.sessions.models import Session

class MetaData(models.Model):
	"""Information about job."""
	
	title = models.CharField(max_length=120, blank=True, null=True)
	description = models.TextField(blank=True, null=True)
	email = models.EmailField(blank=True, null=True)

class Parameters(models.Model):
	"""Specify job parameters."""

	expand_peptides = models.BooleanField(default=False)
	alpha = models.FloatField(default=0.01)
	realign_ragged = models.BooleanField(default=False)
	match_proteases = models.BooleanField(default=False)
	remove_ambiguous = models.BooleanField(default=True)
	multiple_testing_correction = models.BooleanField(default=True)

class DataTypes(models.Model):
	"""Accepted foreground data formats."""

	data_type = models.CharField(max_length=20)

def unique_filename(instance, filename):
	"""Generate directory and UUID name for submitted file"""
	dir_ = type(instance).__name__
	extension = filename.split('.')[-1]
	return '{}/{}.{}'.format(dir_, uuid.uuid4(), extension)

class DataUpload(models.Model):
	"""Job input data from uploaded file."""
	
	file_path = models.FileField('Uploaded File',
								upload_to=unique_filename,
								blank=True,
								null=True)
	raw_seqences = models.TextField(blank=True, null=True)

class Proteome(models.Model):
	"""Accepted proteome selections."""

	file_path = models.FileField('Uploaded File',
									upload_to=unique_filename,
									blank=True,
									null=True)
	raw_seqences = models.TextField(blank=True, null=True)
	organism = models.CharField(max_length=120,
								blank=True,
								null=True)


class BackgroundUpload(models.Model):
	"""Job background frequencies or proteome upload."""
	
	file_path = models.FileField('Uploaded File',
									upload_to=unique_filename,
									null=True)	
	raw_sequences = models.TextField(null=True)
	proteome = models.ForeignKey('Proteome',
									on_delete=models.SET_NULL,
									null=True)
	regular_expression = models.CharField(max_length=120,
											null=True)
	ignore_positions = models.CharField(max_length=120,
										blank=True,
										null=True)

@receiver(models.signals.post_delete, sender=Proteome)
@receiver(models.signals.post_delete, sender=BackgroundUpload)
@receiver(models.signals.post_delete, sender=DataUpload)
def delete_data_file(sender, instance, **kwargs):
	"""
	Delete data file from filesystem after FileUpload object is deleted.
	
	Parameters:
		sender -- Model sending signal.
		instance -- Instance of signal sending class corresponding to
						deleted file.

	Returns:
		None
	"""
	
	if instance.file_path:
		os.remove(instance.file_path.path)

# Should also include a pre-save function to raise error if filename
# already exists.

class Job(models.Model):
	"""Store job details."""

	submitted = models.DateTimeField(auto_now_add=True)
	completed = models.DateTimeField(null=True)
	data_id = models.ForeignKey('DataUpload',
								on_delete=models.CASCADE,
								blank=True,
								null=True)
	background_id = models.ForeignKey('BackgroundUpload',
										on_delete=models.CASCADE,
										blank=True,
										null=True)
	proteome_id = models.ForeignKey('Proteome',
										on_delete=models.CASCADE,
										blank=True,
										null=True)
	parameters_id = models.ForeignKey('Parameters',
										on_delete=models.CASCADE,
										default=1)
	metadata_id = models.ForeignKey('MetaData',
									on_delete=models.DO_NOTHING,
									blank=True,
									null=True)

class Results(models.Model):
	"""Analysis job results."""

	job_id = models.ForeignKey('Job', on_delete=models.CASCADE)