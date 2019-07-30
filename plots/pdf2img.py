import fitz  # from pyMuPdf package
from PIL import Image, ImageChops


def pdf2png(pdf, zoom_ratio=2, rotate=0, out_format='png'):
    doc = fitz.open(pdf)
    for pg in range(doc.pageCount):
        page = doc[pg]
        trans = fitz.Matrix(zoom_ratio, zoom_ratio).preRotate(rotate)
        # create raster image of page (non-transparent)
        pm = page.getPixmap(matrix=trans, alpha=False)
        # write a PNG image of the page
        out = pdf[:-3] + out_format
        pm.writePNG(out)
        return out


def trim_white_around(img):
    path = img
    img = Image.open(img)
    bg = Image.new(img.mode, img.size, img.getpixel((0, 0)))
    diff = ImageChops.difference(img, bg)
    diff = ImageChops.add(diff, diff, 2.0, -90)
    bbox = diff.getbbox()
    if bbox:
        img = img.crop(bbox)
    img.save(path)
    return path


def pdf2img(img, zoom_ratio=2, rotate=0, out_format='png', trim_white=True):
    """
    把第一页的pdf图片抠出来
    :param img:
    :param zoom_ratio:
    :param rotate:
    :param out_format:
    :param trim_white:
    :return:
    """
    img = pdf2png(img, zoom_ratio=zoom_ratio, rotate=rotate, out_format=out_format)
    if trim_white:
        img = trim_white_around(img)
    return img


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pdf2img'])
