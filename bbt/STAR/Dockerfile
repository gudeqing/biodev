FROM alpine
RUN wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz && tar -xzf 2.7.8a.tar.gz
ENV PATH=/STAR-2.7.8a/bin/Linux_x86_64_static/:${PATH}
CMD ["echo", "usage: docker run -it <image_id> STAR --help "]
