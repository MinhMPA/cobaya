"""
.. module:: containers

:Synopsis: Functions and scripts to manage container images
:Author: Jesus Torrado

"""

# Python 2/3 compatibility
from __future__ import absolute_import
from __future__ import division

# Global
import os
import logging
from io import StringIO
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
import uuid
import argparse
from textwrap import dedent

# Local
from cobaya.log import logger_setup, HandledException
from cobaya.input import get_modules, load_input
from cobaya.yaml_custom import yaml_dump
from cobaya.install import install
from cobaya.conventions import _modules_path, _products_path, _code, _data
from cobaya.conventions import _requirements_file, _help_file

logger_setup()
log = logging.getLogger(__name__)

requirements_file_path = os.path.join(_modules_path, _requirements_file)
help_file_path = os.path.join(_modules_path, _help_file)

base_recipe = ur"""
# OS -------------------------------------------------------------------------
FROM ubuntu:xenial
# POST -----------------------------------------------------------------------
RUN sed -i 's/# \(.*multiverse$\)/\1/g' /etc/apt/sources.list && \
    apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y \
      autoconf automake make gcc-6-base \
      libopenblas-base liblapack3 liblapack-dev libcfitsio-dev \
      python python-pip git wget
RUN pip install pip pytest-xdist matplotlib --upgrade
# Prepare environment and tree for modules -----------------------------------
ENV CONTAINED=TRUE
ENV COBAYA_MODULES %s
ENV COBAYA_PRODUCTS %s
RUN mkdir $COBAYA_MODULES && \
    mkdir $COBAYA_PRODUCTS
# COBAYA  --------------------------------------------------------------------
# getdist fork (it will be an automatic requisite in the future)
RUN pip install git+https://github.com/JesusTorrado/getdist/\#egg=getdist
RUN cd $COBAYA_MODULES && git clone https://github.com/JesusTorrado/cobaya.git && \
    cd $COBAYA_MODULES/cobaya && pip install -e .
# Compatibility with singularity ---------------------------------------------
RUN ldconfig
RUN echo "Base image created."
"""%(_modules_path, _products_path)

MPI_recipe = {
    "docker": dedent(u"""
    # MPI -- NERSC: must be MPICH installed in user space
    # http://www.nersc.gov/users/software/using-shifter-and-docker/using-shifter-at-nersc/
    RUN cd /tmp && wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz && \
        tar xvzf mpich-3.2.tar.gz && cd /tmp/mpich-3.2 && ./configure && make -j4 && \
        make install && make clean && rm /tmp/mpich-3.2.tar.gz"""),
    "singularity": u"""apt-get install -y openmpi-bin"""}


def image_help(engine):
    e = engine.lower()
    assert e in ("singularity", "docker"), e + " not valid."
    mount_data = {"docker": "-v [/cluster/path/to/data]:/modules/data:rw",
                  "singularity": "--bind [/cluster/path/to/data]:/modules/data"}
    mount_products = {"docker": "-v [/cluster/path/to/products]:/products:rw",
                      "singularity": "--bind [/cluster/path/to/producs]:/products"}
    mount_tmp = {"docker": "-v /tmp:/products:rw",
                 "singularity": "--bind /tmp:/products"}
    pre_prepare = {"docker": " ".join(["?????", mount_data[e]]),
                   "singularity": " ".join(
                       ["singularity exec", mount_data[e], "[image_file]"])}
    pre_run = {"docker": " ".join(["?????", mount_data[e], mount_products[e]]),
               "singularity": " ".join(
                   [mount_data[e], mount_products[e], "[image_file]"])}
    pre_shell = {"docker": " ".join(["?????", mount_data[e], mount_tmp[e]]),
                 "singularity": " ".join(
                     ["singularity shell", mount_data[e], mount_tmp[e], "[image_file]"])}
    return dedent("""
        This is a %s image for Cobaya.

        To check the modules installed in the container, take a look at '%s'.

        Make sure that you have created a 'data' and a 'products' folder in your cluster.

        To prepare the data needed for the container, while in the cluster, run:

            $ %s cobaya-prepare-data  # --force

        To run a sample with an input file 'somename.yaml', send to you cluster scheduler:

            $ mpi[run|exec] [mpi options] %s cobaya-run somename.yaml

        To open a terminal in the container, for testing purposes do:

            $ %s

        Have fun!
        """%(engine.title(), requirements_file_path, pre_prepare[engine.lower()],
             pre_run[engine.lower()], pre_shell[engine.lower()]))


def get_docker_client():
    try:
        import docker
    except ImportError:
        log.error("The Python Docker interface not installed: do 'pip install docker'.")
        raise HandledException
    return docker.from_env(version="auto")


def create_base_image():
    log.info("Creating base image...")
    dc = get_docker_client()
    with StringIO(base_recipe) as stream:
        dc.images.build(fileobj=stream, tag="cobaya/base:latest", nocache=True)
    log.info("Base image created!")


def create_docker_image(filenames):
    log.info("Creating Docker image...")
    dc = get_docker_client()
    modules = yaml_dump(get_modules(*[load_input(f) for f in filenames])).strip()
    echos_reqs = "RUN "+" && \\ \n    ".join(
        [r'echo "%s" >> %s'%(block, requirements_file_path)
         for block in modules.split("\n")])
    echos_help = "RUN "+" && \\ \n    ".join(
        [r'echo "%s" >> %s'%(line, help_file_path)
         for line in image_help("docker").split("\n")])
    recipe = ur"""
    FROM cobaya/base:latest
    %s
    %s
    RUN cobaya-install %s --path %s --just-code
    %s
    CMD ["cat", "%s"]
    """ % (MPI_recipe["docker"], echos_reqs, requirements_file_path, _modules_path,
           echos_help, help_file_path)
    image_name = "cobaya:"+uuid.uuid4().hex[:6]
    stream = StringIO(recipe)
    dc.images.build(fileobj=stream, tag=image_name)
    stream.close()
    log.info("Docker image '%s' created!", image_name)


def create_singularity_image(*filenames):
    log.info("Creating Singularity image...")
    modules = yaml_dump(get_modules(*[load_input(f) for f in filenames])).strip()
    echos_reqs = "\n".join(['        echo "%s" >> %s'%(block, requirements_file_path)
                            for block in modules.split("\n")])
    recipe = dedent("""
        Bootstrap: docker
        From: cobaya/base:latest
        %%post
          %s
          %s
          cobaya-install %s --path %s --just-code
          mkdir %s
        %%help
        %s
        """ % (MPI_recipe["singularity"], echos_reqs, requirements_file_path,
               _modules_path, os.path.join(_modules_path, _data),
               "\n        ".join(image_help("singularity").split("\n")[1:])))
    print image_help("singularity")
    exit()
    with NamedTemporaryFile(delete=False) as recipe_file:
        recipe_file.write(recipe)
        recipe_file_name = recipe_file.name
    image_name = "cobaya_"+uuid.uuid4().hex[:6]+".simg"
    process_build = Popen(["singularity", "build", image_name, recipe_file_name],
                          stdout=PIPE, stderr=PIPE)
    out, err = process_build.communicate()
    if process_build.returncode:
        log.info(out)
        log.info(err)
        log.error("Image creation failed! See error message above.")
        raise HandledException
    log.info("Singularity image '%s' created!", image_name)


# Command-line scripts ###################################################################

def create_image_script():
    parser = argparse.ArgumentParser(description=(
        "Cobaya's tool for preparing Docker (for Shifter) and Singularity images."))
    parser.add_argument("files", action="store", nargs="+", metavar="input_file.yaml",
                        help="One or more input files.")
    group_type = parser.add_mutually_exclusive_group(required=True)
    group_type.add_argument("-d", "--docker", action="store_const", const="docker",
                            help="Create a Docker image (for Shifter).", dest="type")
    group_type.add_argument("-s", "--singularity", action="store_const", dest="type",
                            const="singularity", help="Create a Singularity image.")
    arguments = parser.parse_args()
    if arguments.type == "docker":
        create_docker_image(arguments.files)
    elif arguments.type == "singularity":
        create_singularity_image(*arguments.files)


def prepare_data_script():
    if "CONTAINED" not in os.environ:
        log.error("This command should only be run within a container. "
                  "Run 'cobaya-install' instead.")
        raise HandledException
    parser = argparse.ArgumentParser(
        description="Cobaya's installation tool for the data needed by a container.")
    parser.add_argument("-f", "--force", action="store_true", default=False,
                        help="Force re-installation of apparently installed modules.")
    arguments = parser.parse_args()
    try:
        info = load_input(requirements_file_path)
    except IOError:
        log.error("Cannot find the requiremets file. This should not be happening.")
        raise HandledException
    install(info, path=_modules_path, force=arguments.force,
            **{_code: False, _data: True})
