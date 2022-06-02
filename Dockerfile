FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/wf-base:fbe8-main

# Its easy to build binaries from source that you can later reference as
# subprocesses within your workflow.
RUN curl -L https://github.com/beliveau-lab/NUPACK/archive/refs/tags/v4.0.0.23.zip -o nupack-4.0.0.23.zip &&\
    unzip nupack-4.0.0.23.zip &&\
    mv nupack-4.0.0.23 /root/nupack

WORKDIR /root/nupack
RUN python3 -m pip install -U nupack -f package


RUN python3 -m pip install -U pip &&\
    python3 -m pip install -U matplotlib jupyterlab


# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
RUN  sed -i 's/latch/wf/g' flytekit.config
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch
WORKDIR /root
