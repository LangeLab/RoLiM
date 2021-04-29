import os
from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse


def index(request):
    return render(request, 'frontend/index.html')


def generate_example_file_response(path):
    file_path = os.path.join(settings.MEDIA_ROOT, path)
    with open(file_path, 'rb') as example_file:
        response = HttpResponse(
            example_file.read(),
            content_type="text/plain"
        )
        response['Content-Disposition'] = (
            'inline; filename=' + os.path.basename(file_path)
        )
    
    return response


def prealigned_text_file(request):
    path = 'examples/prealigned_text_file.txt'
    response = generate_example_file_response(path)

    return response


def text_file_peptide_list(request):
    path = 'examples/text_file_peptide_list.txt'
    response = generate_example_file_response(path)

    return response


def compound_residue_file(request):
    path = 'examples/compound_residue_file_example.txt'
    response = generate_example_file_response(path)

    return response