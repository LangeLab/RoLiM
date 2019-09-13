from rest_framework import serializers
from tool.models import (Job, DataFormat, BackgroundFormat)


class JobSerializer(serializers.ModelSerializer):
	class Meta:
		model = Job
		fields = ('title', 'description', 'email', 'dataupload',
					'backgroundupload', 'alpha', 'submitted', 'completed')
		read_only_fields = ('submitted', 'completed',)


class DataFormatSerializer(serializers.ModelSerializer):
	class Meta:
		model = DataFormat
		fields = '__all__'


class BackgroundFormatSerializer(serializers.ModelSerializer):
	class Meta:
		model = BackgroundFormat
		fields = '__all__'