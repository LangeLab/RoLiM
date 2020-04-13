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
			'position_specific',
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
		extra_kwargs = {
			'position_specific': {'default': True, 'initial': True},
			'extend_sequences': {'default': False, 'initial': False},
			'center_sequences': {'default': True, 'initial': True},
			'multiple_testing_correction': {'default': True, 'initial': True},
			'positional_weighting': {'default': True, 'initial': True},
			'compound_residues': {'default': True, 'initial': True},
			'compound_residue_decomposition': {'default': True, 'initial': True}
		}


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