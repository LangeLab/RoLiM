import os
import tempfile

from django.conf import settings

from tool.models import Job
from logic import sequences, pattern_extraction, merops_connector, vis

DEFAULTS = os.path.join(settings.MEDIA_ROOT, 'defaults')
TEMP = os.path.join(settings.MEDIA_ROOT, 'temp')


def new_job(jobcode):

########################################################################
## INPUT HANDLING LOGIC
########################################################################
    job = Job.objects.get(jobcode=jobcode)

    # Currently only handling FASTA-formatted background.
    if job.backgroundupload:
        # Load context sequences from uploaded FASTA.
        context = sequences.import_fasta(job.backgroundupload.name)
    else:
        # Load context from default Swiss-Prot copy.
        context = sequences.import_fasta(os.path.join(DEFAULTS, 'uniprot.fasta'))
    
    # Generate new Background instance.
    background = sequences.Background(
        context['sequence'].tolist()
    )

    # Figure out file type and run correct pre-processing steps.
    dataformat = job.dataformat_id

    # Create sequence data frame from accepted input file format.
    if dataformat == 1:
        # Load input sequences from prealigned txt file.
        with open(job.dataupload.name) as seqs:
            sequence_df = sequences.sequences_to_df(seqs)
    elif dataformat == 2:
        # Load input sequences from prealigned FASTA file.
        pass
    elif dataformat == 3:
        # Load input sequences from txt peptide list.
        pass
    elif dataformat == 5:
        # Load input sequences from FASTA peptide list.
        pass
    elif dataformat == 6:
        # Load input sequences from text field.
        pass
    elif dataformat == 7:
        # Load input sequences from MaxQuant evidence.txt.
        pass

########################################################################
## NORMAL PATTERN EXTRACTION WORKFLOW
########################################################################
    # Convert sequence data frame to 3D Numpy array.
    sequence_tensor = sequences.vectorize_sequences(sample, background)
    
    kwargs = {
        'max_depth': job.position_limit,
        'p_value_cutoff': job.p_value_cutoff,
        'occurrences': job.minimum_occurrences,
        'fold_change_cutoff': job.fold_change_cutoff,
        'multiple_testing_correction': job.multiple_testing_correction,
        'positional_weighting': job.positional_weighting,
    }

    # Run pattern extraction.
    patterns = logic.PatternContainer(sequence_tensor, background, **kwargs)

########################################################################
## OUTPUT HANDLING LOGIC
########################################################################
# Post-process extracted patterns.
    # Generate logo maps and save to job output directory.
    with tempfile.NamedTemporaryFile() as temp:
        pass
    # Generate sequence txt files and save to job output directory.

    # Add logomap and sequence files to Patterns.
    
# Generate summary outputs.
    # Protease pattern heatmap and save to job output directory.

    # Add heatmap file to Output.

# Mark job as complete.