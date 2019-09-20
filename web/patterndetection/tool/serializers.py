from rest_framework import serializers
from patterndetection.tool.models import (
	Job, ForegroundFormat, ContextFormat, ExtensionDirection
)

class JobSerializer(serializers.ModelSerializer):
	class Meta:
		model = Job
		fields = (
			'email',
			'title',
			'description',
			'foreground_data',
			'foregroundformat',
			'context_data',
			'contextformat',
			'p_value_cutoff',
			'minimum_occurrences',
			'fold_change_cutoff',
			'max_depth',
			'extend_sequences',
			'extension_direction',
			'width',
			'center_sequences',
			'multiple_testing_correction',
			'positional_weighting',
			'compound_residues',
			'compound_residue_decomposition',
		)
		read_only_fields = ('submitted', 'completed',)


class ForegroundFormatSerializer(serializers.ModelSerializer):
	class Meta:
		model = ForegroundFormat
		fields = '__all__'


class ContextFormatSerializer(serializers.ModelSerializer):
	class Meta:
		model = ContextFormat
		fields = '__all__'


class ExtensionDirectionSerializer(serializers.ModelSerializer):
	class Meta:
		model = ExtensionDirection
		fields = '__all__'