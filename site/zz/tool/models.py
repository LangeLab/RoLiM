import os
import uuid

from django.db import models
from django.contrib.sessions.models import Session
from django.dispatch import receiver


def unique_filename(instance, filename):
	"""Generate directory and UUID name for submitted file"""
	dir_ = instance.jobcode
	extension = filename.split('.')[-1]
	return '{}/{}.{}'.format(dir_, uuid.uuid4(), extension)


class Job(models.Model):
	"""Analysis job."""

	DEFAULT_ALPHA = 0.01
	DEFAULT_MINIMUM_OCCURRENCES = 20
	DEFAULT_FC = 1
	DEFAULT_DATAFORMAT = 1
	DEFAULT_BACKGROUNDFORMAT = 1

	jobcode = models.CharField(max_length=32, primary_key=True)
	session = models.ForeignKey(Session,
									on_delete=models.CASCADE)
	submitted = models.DateTimeField(auto_now_add=True)
	completed = models.DateTimeField(blank=True, null=True)
	title = models.CharField(max_length=120,
								blank=True,
								null=True)
	description = models.TextField(blank=True, null=True)
	email = models.EmailField(blank=True, null=True)
	dataupload = models.FileField('Uploaded File',
								upload_to=unique_filename,
								blank=True,
								null=True)
	dataformat = models.ForeignKey('DataFormat',
									on_delete=models.CASCADE,
									default=DEFAULT_DATAFORMAT)
	proteomeupload = models.FileField('Uploaded File',
								upload_to=unique_filename,
								blank=True,
								null=True)
	backgroundupload = models.FileField('Uploaded File',
								upload_to=unique_filename,
								blank=True,
								null=True)
	backgroundformat = models.ForeignKey('BackgroundFormat',
										 on_delete=models.CASCADE,
										 default=DEFAULT_BACKGROUNDFORMAT)
	alpha = models.FloatField(default=DEFAULT_ALPHA)
	multiple_testing_correction = models.BooleanField(default=True)
	minimum_occurrences = models.IntegerField(default=DEFAULT_MINIMUM_OCCURRENCES)
	fold_change_cutoff = models.FloatField(default=DEFAULT_FC)
	expand_peptides = models.BooleanField(default=False)
	center_positions = models.BooleanField(default=True)
	compound_residues = models.BooleanField(default=True)
	position_limit = models.IntegerField(blank=True, null=True)
	positional_weighting = models.BooleanField(default=True)
	ignore_background_positions = models.CharField(max_length=120,
											blank=True,
											null=True)
	background_regular_expression = models.CharField(max_length=120,
														blank=True,
														null=True)


class Pattern(models.Model):
	"""Pattern details."""
	job = models.ForeignKey('Job', on_delete=models.CASCADE)
	logomap = models.FileField(upload_to=unique_filename)
	sequences = models.FileField(upload_to=unique_filename)
	pattern_string = models.CharField(max_length=120)
	order = models.IntegerField()


class Figure(models.Model):
	"""Default output figures."""
	job = models.ForeignKey('Job', on_delete=models.CASCADE)
	figure = models.FileField(upload_to=unique_filename)
	data = models.FileField(upload_to=unique_filename)
	figuretype = models.ForeignKey('FigureType', on_delete=models.CASCADE)


class FigureType(models.Model):
	"""Descriptions of default output figures."""
	description = models.TextField()


class BackgroundFormat(models.Model):
	"""Accepted background sequence set types."""

	background_format = models.CharField(max_length=120)


class DataFormat(models.Model):
	"""Accepted input data formats."""

	data_format = models.CharField(max_length=120)


@receiver(models.signals.post_delete, sender=Job)
def delete_data_file(sender, instance, **kwargs):
	"""
	Delete data file from filesystem after file upload object is deleted.
	
	Parameters:
		sender -- Model sending signal.
		instance -- Instance of signal sending class corresponding to
						deleted file.

	Returns:
		None
	"""
	
	data_files = [
		instance.dataupload,
		instance.backgroundupload,
		instance.proteomeupload,
	]

	for data_file in data_files:
		if data_file:
			os.remove(data_file)