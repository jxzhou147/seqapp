from django.shortcuts import render

from .functions import sequence_len


def index(request):
    return render(request, 'try/index.html')


def results(request):
    sequence_input = request.POST['sequence']
    sequence_length = sequence_len(sequence_input)
    context = {'output': sequence_length}
    return render(request, 'try/results.html', context)
