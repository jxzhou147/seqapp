from django.shortcuts import render

from .reproduction import run


def index(request):
    return render(request, 'try/index.html')


def results(request):
    sequence_input = request.POST['sequence']
    NET_CHARGE, CDH3_HI, FabNetCharge, FvCSP, HISum = run(sequence_input)
    context = {'net_charge': NET_CHARGE,
            'HIndex': CDH3_HI,
            'FabNetCharge': FabNetCharge,
            'FvCSP': FvCSP,
            'HISum': HISum}
    return render(request, 'try/results.html', context)
