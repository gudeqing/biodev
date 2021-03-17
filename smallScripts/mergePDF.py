import glob
import sys
import os
from PyPDF2 import PdfFileWriter, PdfFileReader, PdfFileMerger


def merger(files:list, out):
    pdf_merger = PdfFileMerger()

    for ind, path in enumerate(files):
        pdf_merger.append(path)
        title = os.path.basename(path).split('.', 1)[0]
        pdf_merger.addBookmark(title, ind, parent=None)

    pdf_merger.setPageLayout(layout='/TwoColumnLeft')

    with open(out, 'wb') as fileobj:
        pdf_merger.write(fileobj)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['merger'])

