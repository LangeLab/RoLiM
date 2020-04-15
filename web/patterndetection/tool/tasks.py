import os
import re
import shutil
import smtplib

from django.conf import settings
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

from patterndetection.tool.models import Job
from zigzag import sequences, pattern_extraction

DEFAULTS = os.path.join(settings.MEDIA_ROOT, 'defaults')
TEMP = os.path.join(settings.MEDIA_ROOT, 'temp')


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
    compound_residues = job.compound_residues
    compound_residue_decomposition = job.compound_residue_decomposition
    position_specific = job.position_specific

    # Set terminal based on sequence extension direction.
    if extension_direction == 1:
        terminal = 'n'
    elif extension_direction == 2:
        terminal = 'c'
    elif extension_direction == 3:
        terminal = 'both'
    else:
        terminal = None

    foreground_file_name = os.path.join(settings.MEDIA_ROOT, foreground_data.name)
    
    # Set email login parameters.
    username = 'tsmithdmr@gmail.com'
    with open('gmail_app_password') as app_password:
        password = app_password.read().strip()
    
    # Generate new email message and specify details.
    msg = MIMEMultipart('alternative')
    msg['Subject'] = 'RoLiM Results: {}'.format(title)
    msg['From'] = username
    msg['To'] = email

    try:
        # Currently only handling FASTA-formatted background.
        if context_data:
            # Load context sequences from uploaded FASTA.
            context = sequences.import_fasta(
                os.path.join(settings.MEDIA_ROOT, context_data.name)
            )
            # Generate new Background instance.
            if compound_residues:
                background = sequences.Background(
                    context['sequence'].tolist(),
                    width=width,
                    position_specific=position_specific
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
            context = sequences.import_fasta(os.path.join(DEFAULTS, 'uniprot.fasta'))
            # Generate new Background instance.
            if compound_residues:
                background = sequences.Background(
                    context['sequence'].tolist(),
                    width=width,
                    position_specific=position_specific,
                    precomputed=os.path.join(DEFAULTS, 'swissprot_human_background_{}.csv'.format(width))
                )
            else:
                background = sequences.Background(
                    context['sequence'].tolist(),
                    width=width,
                    position_specific=position_specific,
                    compound_residues=None,
                    precomputed=os.path.join(DEFAULTS, 'swissprot_human_background.csv_{}'.format(width))
                )

        # Load sequences from Job submission data set.
        if foreground_format == 1:
            # Load input sequences from prealigned txt file.
            sample = sequences.load_prealigned_file(
                foreground_file_name,
                background,
                center=center_sequences
            )
        elif foreground_format == 2:
            # Load input sequences from prealigned FASTA file.
            sample = sequences.load_prealigned_fasta(
                foreground_file_name,
                background,
                center=center_sequences
            )
        elif foreground_format == 3:
            # Load input sequences from txt peptide list.
            sample = sequences.load_peptide_list_file(
                foreground_file_name,
                context,
                background,
                center=center_sequences,
                width=width,
                terminal=terminal
            )
        elif foreground_format == 5:
            # Load input sequences from FASTA peptide list.
            sample = sequences.load_fasta_peptides(
                foreground_file_name,
                context,
                background,
                center=center_sequences,
                width=width,
                terminal=terminal
            )
        elif foreground_format == 6:
            # Load input sequences from text field.
            sample = sequences.load_prealigned_field(
                foreground_file_name,
                background,
                center=center_sequences
            )
        elif foreground_format == 7:
            # Load input sequences from MaxQuant evidence.txt.
            sample = sequences.load_maxquant_evidence_file(
                foreground_file_name,
                context,
                background,
                width=width
            )

        output_directory = os.path.join(settings.MEDIA_ROOT, jobcode, 'results')

        # Generate log file.
        log_file_path = os.path.join(output_directory, 'log.txt')
        with open(log_file_path, 'w') as log_file:
            log_file.write('Title:  {}\n'.format(title))
            log_file.write('Description:  {}\n'.format(description))
            log_file.write('Email:  {}\n'.format(email))
            log_file.write('Submitted:  {}\n'.format(submitted))
            log_file.write('Foreground file:  {}\n'.format(foreground_filename))
            log_file.write(
                'Foreground format:  {}\n'.format(
                    ForegroundFormat.get(id=foreground_format)
                )
            )
            log_file.write(
                'Context file:  {}\n'.format(
                    context_filename if context_filename
                    else "Swiss-Prot human proteome"
                )
            )
            log_file.write(
                'Context format:  {}\n'.format(
                    ContextFormat.get(id=context_format)
                )
            )
            log_file.write('P-value cutoff:  {}\n'.format(p_value_cutoff))
            log_file.write('Minimum occurrences:  {}\n'.format(minimum_occurrences))
            log_file.write('Fold change cutoff:  {}\n'.format(fold_change_cutoff))
            log_file.write('Max depth:  {}\n'.format(max_depth))
            log_file.write('Sequence extension:  {}\n'.format(extend_sequences))
            log_file.write('Extension direction:  {}\n'.format(extension_direction))
            log_file.write('Width:  {}\n'.format(width))
            log_file.write('Centered sequences:  {}\n'.format(center_sequences))
            log_file.write(
                'Multiple testing correction:  {}\n'.format(multiple_testing_correction)
            )
            log_file.write(
                'Positional weighting:  {}\n'.format(positional_weighting)
            )
            log_file.write('Compound residue detection:  {}\n'.format(compound_residues))
            log_file.write(
                'Compound residue decomposition:  {}\n'.format(compound_residue_decomposition)
            )         

        # Run pattern extraction analysis and generate outputs.
        patterns = pattern_extraction.PatternContainer(
            sample,
            background,
            title,
            output_directory,
            max_depth=max_depth,
            p_value_cutoff=p_value_cutoff,
            minimum_occurrences=minimum_occurrences,
            fold_change_cutoff=fold_change_cutoff,
            multiple_testing_correction=multiple_testing_correction,
            positional_weighting=positional_weighting,
            allow_compound_residue_decomposition=compound_residue_decomposition
        )
        if width % 2 == 0:
            patterns.post_processing()
        else:
            patterns.post_processing(proteolysis_data=False)

        # Compress outputs as zip archive.
        shutil.make_archive(
            output_directory,
            'zip',
            os.path.join(settings.MEDIA_ROOT, jobcode),
            'results'
        )

        # Prepare email.
        html = 'Please find attached, the results of your pattern detection analysis.'
        msg_body = MIMEText(html, 'html')
        msg.attach(msg_body)
        
        # Attach results archive to email message.
        filename = 'results.zip'
        zip_path = os.path.join(settings.MEDIA_ROOT, jobcode, 'results.zip')
        with open(zip_path, 'rb') as attachment:
            p = MIMEBase('application', 'octet-stream')
            p.set_payload((attachment).read())
            encoders.encode_base64(p)
            p.add_header('Content-Disposition', "attachment; filename= %s" % filename)
            msg.attach(p)
    except Exception as e:
        html = (
            'Something went wrong with your analysis:<br /><br />{}'.format(re.escape(e))
            + '<br /><br />Please check that the options you selected match the'
            + ' format of your data, and that the format of your data matches a'
            + ' format supported by our tool.<br /><br />Thank you!'
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
    job.delete()