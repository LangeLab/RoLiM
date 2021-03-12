from rest_framework import serializers
from patterndetection.tool.models import (
	Job, ForegroundFormat, ContextFormat, ExtensionDirection,
	RedundancyLevel, OriginalRowMerge,
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
			'foreground_filename',
			'contextformat',
			'context_data',
			'context_filename',
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
			'compoundresidue_file',
			'compoundresidue_filename',
			'require_context_id',
			'redundancylevel_id',
			'originalrowmerge_id',
			'first_protein_only',
			'cluster_sequences',
		)
		read_only_fields = ('submitted', 'completed',)
		extra_kwargs = {
			'position_specific': {'default': True, 'initial': True},
			'extend_sequences': {'default': False, 'initial': False},
			'center_sequences': {'default': True, 'initial': True},
			'multiple_testing_correction': {'default': True, 'initial': True},
			'positional_weighting': {'default': True, 'initial': True},
			'compound_residues': {'default': True, 'initial': True},
			'compound_residue_decomposition': {'default': True, 'initial': True},
			'require_context_id': {'default': True, 'initial': True},
			'cluster_sequences': {'default': True, 'initial': True},
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

class RedundancyLevelSerializer(serializers.ModelSerializer):
	class Meta:
		model = RedundancyLevel
		fields = '__all__'


class OriginalRowMergeSerializer(serializers.ModelSerializer):
	class Meta:
		model = OriginalRowMerge
		fields = '__all__'