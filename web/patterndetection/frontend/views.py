from django.shortcuts import render

def index(request):
	return render(request, 'frontend/index.html')


def splash(request):
	return render(request, 'frontend/splash.html')
