import os
from os.path import join, exists, isdir, dirname
from commands import getstatusoutput as getoutput

def rsystem(command):
    """
    Executes command, raise an error if it fails and send status to logfile
    :param command: command to be executed
    :return: 
    """
    print("\nEXE:", command)
    
    if os.system(command) != 0:
        print("Error executing:", command)
        raise TypeError(command)
    else:
        print("Done executing:", command)

def download_ehm():

    source = 'https://github.com/txusser/ExactHistogramSpecification.git'

    icom = 'git clone %s exact_hm' % source
    rsystem(icom)

def unpack_example():

    command = 'unzip resources.zip'
    rsystem(command)

def install_soap():
    """
    Execute installation of the required packages and libraries
    :return: 
    """
    print("Insatalling SOAP")

    icom = 'sudo apt update'
    rsystem(icom)

    # Install and upgrade PIP
    icom = 'sudo apt -y install python-pip'
    rsystem(icom)

    icom = 'sudo apt -y install unzip'
    rsystem(icom)

    icom = 'pip install --upgrade pip'
    rsystem(icom)

    # Install numpy
    icom = 'sudo  pip install numpy'
    rsystem(icom)

    # Install Nibabel
    icom = 'sudo pip install nibabel'
    rsystem(icom)

    icom = 'sudo pip install scipy'
    rsystem(icom)


install_soap()
download_ehm()
unpack_example()
