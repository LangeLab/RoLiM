from rest_framework import serializers
from tool.models import DataUpload


class DataUploadSerializer(serializers.HyperlinkedModelSerializer):
	"""Serialize uploaded file."""

	class Meta:
		model = DataUpload
		read_only_fields = ('file_path')

class BackgroundUploadSerializer(serializers.HyperlinkedModelSerializer):
	"""Serialize uploaded file."""

	class Meta:
		model = BackgroundUpload
		read_only_fields = ('file_path')

class DataTypesSerializer(serializers.HyperlinkedModelSerializer):
	pass

class ParametersSerializer(serializers.HyperlinkedModelSerializer):
	pass

class MetaDataSerializer(serializers.HyperlinkedModelSerializer):
	pass

class ProteomeSerializer(serializers.HyperlinkedModelSerializer):
	pass

class JobSerializer(serializers.HyperlinkedModelSerializer):
	pass
