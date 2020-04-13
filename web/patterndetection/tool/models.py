import os
import shutil
import uuid

from django.conf import settings
from django.db import models
from django.contrib.sessions.models import Session
from django.dispatch import receiver


def unique_filename(instance, filename):
	"""Generate directory and UUID name for submitted file"""
	dir_ = instance.jobcode
	extension = filename.split('.')[-1]
	return '{}/{}.{}'.format(dir_, uuid.uuid4(), extension)


class Job(models.Model):
	"""Analysis input data and metadata."""

	DEFAULT_P_VALUE_CUTOFF = 0.001
	DEFAULT_MIN_OCCURRENCES = 2
	DEFAULT_FC = 1.0
	DEFAULT_DATAFORMAT = 1
	DEFAULT_CONTEXTFORMAT = 1
	DEFAULT_EXTENSION_DIRECTION = 1

	jobcode = models.CharField(max_length=32, primary_key=True)
	session = models.ForeignKey(Session, on_delete=models.CASCADE)
	submitted = models.DateTimeField(auto_now_add=True)
	completed = models.DateTimeField(blank=True, null=True)
	email = models.EmailField()
	title = models.CharField(max_length=120)
	description = models.TextField(blank=True, null=True)
	foreground_data = models.FileField('Uploaded File',
										upload_to=unique_filename,
										blank=True,
										null=True)
	foregroundformat = models.ForeignKey('ForegroundFormat',
											on_delete=models.CASCADE,
											default=DEFAULT_DATAFORMAT)
	context_data = models.FileField('Uploaded File',
										upload_to=unique_filename,
										blank=True,
										null=True)
	contextformat = models.ForeignKey('ContextFormat',
										on_delete=models.CASCADE,
										default=DEFAULT_CONTEXTFORMAT)
	p_value_cutoff = models.FloatField(default=DEFAULT_P_VALUE_CUTOFF)
	minimum_occurrences = models.IntegerField(blank=True,
												default=DEFAULT_MIN_OCCURRENCES)
	fold_change_cutoff = models.FloatField(blank=True, default=DEFAULT_FC)
	max_depth = models.IntegerField(blank=True, null=True)
	extend_sequences = models.BooleanField(blank=True, default=False)
	extension_direction = models.ForeignKey('ExtensionDirection',
												on_delete=models.CASCADE,
												default=DEFAULT_EXTENSION_DIRECTION)
	width = models.IntegerField(blank=True, default=8)
	center_sequences = models.BooleanField(blank=True, default=True)
	multiple_testing_correction = models.BooleanField(blank=True, default=True)
	positional_weighting = models.BooleanField(blank=True, default=True)
	compound_residues = models.BooleanField(blank=True, default=True)
	compound_residue_decomposition = models.BooleanField(blank=True, default=True)
	position_specific = models.BooleanField(blank=True, default=True)


class ForegroundFormat(models.Model):
	"""Supported foreground data set format options."""
	foreground_format = models.CharField(max_length=120)


class ContextFormat(models.Model):
	"""Supported context data set format options."""
	context_format = models.CharField(max_length=120)


class ExtensionDirection(models.Model):
	"""Supported sequence extension directions."""
	direction = models.CharField(max_length=120)


@receiver(models.signals.post_delete, sender=Job)
def delete_job_data(sender, instance, **kwargs):
	"""
	Delete data file from filesystem after file upload object is deleted.
	
	Parameters:
		sender -- Model sending signal.
		instance -- Instance of signal sending class corresponding to
						deleted file.

	Returns:
		None
	"""
	
	job_directory = os.path.join(settings.MEDIA_ROOT, instance.jobcode)
	try:
		shutil.rmtree(job_directory)
	except:
		print('Error while removing job files.')