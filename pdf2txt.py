import PyPDF2


def pdf2txt(pdf, page_range=(0, 1), out='pdf.txt'):
    with open(pdf, 'rb') as f, open(out, 'w') as f2:
        pdfReader = PyPDF2.PdfFileReader(f)
        print(" No. Of Pages :", pdfReader.numPages)
        for i in range(page_range[0], page_range[1]):
            pageObject = pdfReader.getPage(i)
            f2.write(pageObject.extractText())


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pdf2txt'])

