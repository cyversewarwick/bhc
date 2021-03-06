#run as:
#	docker build -t bhc .

#grab existing R image
FROM r-base:latest

#need Python for R's argparse port and subsequent BiNGO/MEME postprocessing
RUN apt-get -y update && apt-get -y upgrade
RUN apt-get -y install python3 python3-numpy ttf-bitstream-vera

#copy over scripts
RUN mkdir /scripts
COPY scripts /scripts

#R-side setup
RUN Rscript /scripts/setup.R

#set up analysis crash text file
RUN apt-get -y install git
RUN git clone https://github.com/cyversewarwick/analysis_crash.git

#that's it, ready to use
ENTRYPOINT ["bash", "/scripts/bhc_tarwrapper.sh"]
