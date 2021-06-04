FROM alpine
RUN wget http://opengene.org/fastp/fastp && chmod a+x ./fastp \
	&& mv fastp /bin/ \
	&& apk add --no-cache bash
CMD ["echo", "usage: docker run -it <image_id> fastp ./fastp -i read1 -I read2 "]
