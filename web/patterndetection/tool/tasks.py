import os
import traceback
import re
import shutil
import smtplib

import pandas as pd

from django.conf import settings
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

from patterndetection.tool.models import (
    Job, ForegroundFormat, ContextFormat, RedundancyLevel, OriginalRowMerge
)
from zigzag import sequences, pattern_extraction

DEFAULTS = os.path.join(settings.MEDIA_ROOT, 'defaults')
TEMP = os.path.join(settings.MEDIA_ROOT, 'temp')

MAX_MEMORY_ALLOWANCE = 24.0


class SpeciesError(Exception): pass


class ContextError(Exception): pass


class ForegroundError(Exception): pass


class ParameterError(Exception): pass


class MemoryLimitExceededError(Exception): pass


def check_memory_limit(context, width):
    """
    Checks if job is likely to cause a memory error. Returns True if
        esimtated peak memory load for an analysis is below max memory
        allowance and False if it exceeds this limit.

    Parameters:
        context -- Pandas DataFrame.
        width -- Int.

    Returns:
        below_memory_limit -- Boolean.
    """

    if sequences.estimate_peak_memory_load(context, width) < MAX_MEMORY_ALLOWANCE:
        below_memory_limit = True
    else:
        below_memory_limit = False
    
    return below_memory_limit


def load_compound_residues(compound_residue_file):
    if compound_residue_file:
        compound_residue_path = os.path.join(
            settings.MEDIA_ROOT,
            compound_residue_file.name
        )
        compound_residues = sequences.load_compound_residue_file(
            compound_residue_path
        )
    else:
        compound_residues = sequences.COMPOUND_RESIDUES

    return compound_residues


def new_job(jobcode):
    """
    Process and run job from user submission.

    Parameters:
        jobcode -- UUID. Primary key of submission in Job table.
    
    Returns:
        None
    """

    ########################################################################
    ## INPUT HANDLING LOGIC
    ########################################################################
    job = Job.objects.get(jobcode=jobcode)

    # Get job data from Job submission row.
    foreground_data = job.foreground_data
    context_data = job.context_data
    compound_residue_file = job.compoundresidue_file

    # Get job metadata from Job submission row.
    email = job.email
    title = job.title
    description = job.description
    submitted = job.submitted
    foreground_format = job.foregroundformat_id
    foreground_filename = job.foreground_filename
    context_format = job.contextformat_id
    context_filename = job.context_filename
    p_value_cutoff = job.p_value_cutoff
    minimum_occurrences = job.minimum_occurrences
    fold_change_cutoff = job.fold_change_cutoff
    max_depth = job.max_depth
    extend_sequences = True if foreground_format not in [1, 2, 6] else False
    extension_direction = job.extension_direction_id
    width = job.width
    center_sequences = job.center_sequences
    multiple_testing_correction = job.multiple_testing_correction
    positional_weighting = job.positional_weighting
    enable_compound_residues = job.compound_residues
    compound_residue_decomposition = job.compound_residue_decomposition
    compound_residue_filename = job.compoundresidue_filename
    position_specific = job.position_specific
    require_context_id = job.require_context_id
    first_protein_only = job.first_protein_only
    redundancy_level = RedundancyLevel.objects.get(id=job.redundancylevel_id).redundancy_level
    original_row_merge = OriginalRowMerge.objects.get(id=job.originalrowmerge_id).original_row_merge
    cluster_sequences = job.cluster_sequences
    
    # Set terminal based on sequence extension direction.
    if extension_direction == 1:
        terminal = 'n'
    elif extension_direction == 2:
        terminal = 'c'
    elif extension_direction == 3:
        terminal = 'both'
    else:
        terminal = None

    foreground_path = os.path.join(settings.MEDIA_ROOT, foreground_data.name)
    
    # Set email login parameters.
    username = 'lange.lab.ubc@gmail.com'
    with open('gmail_app_password') as app_password:
        password = app_password.read().strip()
    
    # Generate new email message and specify details.
    msg = MIMEMultipart('mixed')
    msg['Subject'] = 'RoLiM Results: {}'.format(title)
    msg['From'] = username
    msg['To'] = email

    output_title = re.sub(r'\W+', ' ', title).strip().replace(" ", "_")
    output_directory = os.path.join(settings.MEDIA_ROOT, jobcode, output_title)
    try:
        os.makedirs(output_directory)
    except:
        pass

    # Generate log file.
    log_file_directory = os.path.join(output_directory, 'summary')
    try:
        os.makedirs(log_file_directory)
    except:
        pass
    log_file_path = os.path.join(log_file_directory, 'log.txt')
    with open(log_file_path, 'w') as log_file:
        log_file.write('Title:  {}\n'.format(title))
        log_file.write('Description:  {}\n'.format(description))
        log_file.write('Email:  {}\n'.format(email))
        log_file.write('Submitted:  {}\n'.format(submitted))
        log_file.write('Foreground file:  {}\n'.format(foreground_filename))
        log_file.write(
            'Foreground format:  {}\n'.format(
                ForegroundFormat.objects.get(id=foreground_format).foreground_format
            )
        )
        log_file.write(
            'Context file:  {}\n'.format(
                context_filename if context_filename
                else "Swiss-Prot {} proteome".format(
                    'human' if context_format == 2 else 'mouse'
                )
            )
        )
        log_file.write(
            'Context format:  {}\n'.format(
                ContextFormat.objects.get(id=context_format).context_format
            )
        )
        log_file.write('P-value cutoff:  {}\n'.format(p_value_cutoff))
        log_file.write('Minimum occurrences:  {}\n'.format(minimum_occurrences))
        log_file.write('Fold difference cutoff:  {}\n'.format(fold_change_cutoff))
        log_file.write('Max depth:  {}\n'.format(max_depth))
        log_file.write('Sequence extension:  {}\n'.format(extend_sequences))
        log_file.write('Extension direction:  {}\n'.format(terminal))
        log_file.write('Require protein identifiers:  {}\n'.format(require_context_id))
        log_file.write('Width:  {}\n'.format(width))
        log_file.write('Centered sequences:  {}\n'.format(center_sequences))
        log_file.write('Redundancy elimination level:  {}\n'.format(redundancy_level))
        log_file.write(
            'Merge multiple aligned peptides from same original row:  {}\n'.format(
                original_row_merge
            )
        )
        log_file.write(
            'Multiple testing correction:  {}\n'.format(multiple_testing_correction)
        )
        log_file.write(
            'Positional weighting:  {}\n'.format(positional_weighting)
        )
        log_file.write('Position specific background:  {}\n'.format(position_specific))
        log_file.write('Compound residue detection:  {}\n'.format(enable_compound_residues))
        log_file.write(
            'Compound residue decomposition:  {}\n'.format(compound_residue_decomposition)
        )
        log_file.write(
            'Compound residue file:  {}\n'.format(
                compound_residue_filename if compound_residue_filename else ""
            )
        )

    try:
        try:
            # Currently only handling FASTA-formatted background.
            if context_data:
                # Load context sequences from uploaded FASTA.
                context = sequences.import_fasta(
                    os.path.join(settings.MEDIA_ROOT, context_data.name)
                )

                if not check_memory_limit(context, width):
                    raise MemoryLimitExceededError()

                # Generate new Background instance.
                if enable_compound_residues:
                    compound_residues = load_compound_residues(compound_residue_file)
                    background = sequences.Background(
                        context['sequence'].tolist(),
                        width=width,
                        position_specific=position_specific,
                        compound_residues=compound_residues
                    )
                else:
                    background = sequences.Background(
                        context['sequence'].tolist(),
                        width=width,
                        position_specific=position_specific,
                        compound_residues=None
                    )
            else:
                # Load context from default Swiss-Prot copy.
                if context_format == 2:
                    context = sequences.import_fasta(os.path.join(DEFAULTS, 'uniprot.fasta'))
                elif context_format == 3:
                    context = sequences.import_fasta(os.path.join(DEFAULTS, 'uniprot_Mus_musculus.fasta'))

                if not check_memory_limit(context, width):
                    raise MemoryLimitExceededError()

                # Generate new Background instance.
                try:
                    if context_format == 2:
                        if enable_compound_residues:
                            compound_residues = load_compound_residues(compound_residue_file)
                            background = sequences.Background(
                                context['sequence'].tolist(),
                                width=width,
                                position_specific=position_specific,
                                precomputed=os.path.join(DEFAULTS, 'swissprot_human_background_{}.csv'.format(width)),
                                compound_residues=compound_residues
                            )
                        else:
                            background = sequences.Background(
                                context['sequence'].tolist(),
                                width=width,
                                position_specific=position_specific,
                                compound_residues=None,
                                precomputed=os.path.join(DEFAULTS, 'swissprot_human_background_{}.csv'.format(width))
                        )
                    else:
                        raise SpeciesError()
                except:
                    if enable_compound_residues:
                        compound_residues = load_compound_residues(compound_residue_file)
                        background = sequences.Background(
                            context['sequence'].tolist(),
                            width=width,
                            position_specific=position_specific,
                            compound_residues=compound_residues
                        )
                    else:
                        background = sequences.Background(
                            context['sequence'].tolist(),
                            width=width,
                            position_specific=position_specific,
                            compound_residues=None
                    )
        except Exception as e:
            raise ContextError() from e

        try:
            # Load sequences from Job submission data set.
            if foreground_format == 1:
                # Load input sequences from prealigned txt file.
                samples = sequences.load_prealigned_file(
                    foreground_path,
                    background,
                    center=center_sequences,
                    redundancy_level=redundancy_level,
                    title=title
                )
            ### NOT ENABLED ###
            elif foreground_format == 2:
                # Load input sequences from prealigned FASTA file.
                samples = sequences.load_prealigned_fasta(
                    foreground_path,
                    background,
                    center=center_sequences,
                    redundancy_level=redundancy_level
                )
            elif foreground_format == 3:
                # Load input sequences from txt peptide list.
                samples = sequences.load_peptide_list_file(
                    foreground_path,
                    context,
                    background,
                    center=center_sequences,
                    width=width,
                    terminal=terminal,
                    require_context_id=require_context_id,
                    redundancy_level=redundancy_level,
                    first_protein_only=first_protein_only,
                    original_row_merge=original_row_merge,
                    title=title
                )
            ### NOT ENABLED ###
            elif foreground_format == 5:
                # Load input sequences from FASTA peptide list.
                samples = sequences.load_fasta_peptides(
                    foreground_path,
                    context,
                    background,
                    center=center_sequences,
                    width=width,
                    terminal=terminal,
                    require_context_id=require_context_id,
                    redundancy_level=redundancy_level,
                    first_protein_only=first_protein_only,
                    original_row_merge=original_row_merge
                )
            ### NOT ENABLED ###
            elif foreground_format == 6:
                # Load input sequences from text field.
                samples = sequences.load_prealigned_field(
                    foreground_path,
                    background,
                    center=center_sequences,
                    redundancy_level=redundancy_level
                )
            ### NOT ENABLED ###
            elif foreground_format == 7:
                # Load input sequences from MaxQuant evidence.txt.
                samples = sequences.load_maxquant_evidence_file(
                    foreground_path,
                    context,
                    background,
                    width=width
                )
        except Exception as e:
            raise ForegroundError() from e

        try:
            sample_output_paths = []
            if len(samples.keys()) > 1:
                for sample_name in samples.keys():
                    sample_directory = re.sub(r'\W+', ' ', sample_name).strip().replace(" ", "_")
                    sample_output_path = os.path.join(output_directory, sample_directory)
                    sample_output_paths.append((sample_name, sample_output_path))
            else:
                sample_output_paths.append((list(samples.keys())[0], output_directory))
            
            # Run pattern extraction analysis and generate outputs.
            all_pattern_containers = []
            summary_tables = []
            for sample_name, sample_output_path in sample_output_paths:
                patterns = pattern_extraction.PatternContainer(
                    samples[sample_name],
                    background,
                    str(sample_name),
                    sample_output_path,
                    max_depth=max_depth,
                    p_value_cutoff=p_value_cutoff,
                    minimum_occurrences=minimum_occurrences,
                    fold_change_cutoff=fold_change_cutoff,
                    multiple_testing_correction=multiple_testing_correction,
                    positional_weighting=positional_weighting,
                    allow_compound_residue_decomposition=compound_residue_decomposition
                )
                if (width == 8) and center_sequences:
                    summary_tables.append(
                        patterns.post_processing(
                            cluster_sequences=cluster_sequences,
                            logo_maps=False
                        )
                    )
                else:
                    summary_tables.append(
                        patterns.post_processing(
                            proteolysis_data=False,
                            cluster_sequences=cluster_sequences,
                            logo_maps=False
                        )
                    )
                all_pattern_containers.append(patterns)

            # Add sequence summary table to summary output directory.
            if len(summary_tables) > 1:
                unique_pattern_strings = []
                unique_pattern_list = []
                for pattern_container in all_pattern_containers:
                    for pattern in pattern_container.pattern_list:
                        pattern_string = ''.join(pattern.character_pattern())
                        if pattern_string not in unique_pattern_strings:
                            unique_pattern_strings.append(pattern_string)
                            unique_pattern_list.append(pattern)
                all_original_sequences = pd.concat(
                    [sample.original_sequences for sample in samples.values()]
                )
                summary_table = pattern_extraction.PatternContainer.generate_summary_table(
                    all_original_sequences,
                    unique_pattern_list
                )
                summary_table.to_csv(
                    os.path.join(
                        log_file_directory,
                        '{}_sequence_summary_table.txt'.format(output_title)
                    ),
                    sep='\t',
                    header=True,
                    index=True
                )

            # Compress outputs as zip archive.
            archive_name = re.sub(r'\W+', ' ', title).strip().replace(" ", "_")
            zip_path = os.path.join(settings.MEDIA_ROOT, jobcode, archive_name)
            shutil.make_archive(
                zip_path,
                'zip',
                os.path.join(settings.MEDIA_ROOT, jobcode),
                output_title
            )
            
            # Attach results archive to email message.
            zip_filename = zip_path + '.zip'
            zip_basename = os.path.basename(zip_filename)
            with open(zip_filename, 'rb') as attachment:
                p = MIMEBase('application', 'octet-stream')
                p.set_payload((attachment).read())
                encoders.encode_base64(p)
                p.add_header('Content-Disposition', "attachment; filename= %s" % zip_basename)
                msg.attach(p)

            # Prepare email.
            html = 'Please find attached, the results of your pattern detection analysis.\r\n'
            msg_body = MIMEText(html, 'html')
            msg.attach(msg_body)
            
        except Exception as e:
            raise ParameterError() from e
    except Exception as e:
        exception_type = type(e).__name__
        if exception_type == 'ContextError':
            error_message = (
                'Please check the format and contents of the context'
                + ' file you selected'
            )
        elif exception_type == 'ForegroundError':
            error_message = (
                'Please check the format and contents of the foreground'
                + ' file you provided'
            )
        elif exception_type == 'ParameterError':
            error_message = (
                'Please ensure that the options you selected are'
                + ' consistent with your data'
            )
        elif exception_type == 'MemoryLimitExceededError':
            error_message = (
                'The combination of window width and context data set'
                + ' you selected results in a background which exceeds'
                + '  the total memory allowance for a single analysis.'
                + ' Please reduce the size of one or both in order to'
                + ' try again'
            )
        else:
            error_message = 'review your submission and try again'

        with open(log_file_path, 'rb') as attachment:
            p = MIMEBase('application', 'octet-stream')
            p.set_payload((attachment).read())
            encoders.encode_base64(p)
            p.add_header(
                'Content-Disposition',
                "attachment; filename= %s" % os.path.basename(log_file_path)
            )
            msg.attach(p)

        html = (
            'Something went wrong with your analysis. {}:<br /><br />{}'.format(
                error_message,
                traceback.format_exc()
            )
            + '<br /><br />Thank you!\r\n'
        )
        msg_body = MIMEText(html, 'html')
        msg.attach(msg_body)

    # Start SMTP session using TLS security and login to Gmail
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.connect('smtp.gmail.com', 587)
    server.ehlo()
    server.starttls()
    server.login(username, password)
    server.sendmail(username, email, msg.as_string())
    server.quit()

    # Remove completed job from Job table.
    # job.delete()