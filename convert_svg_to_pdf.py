from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM
from pdf2image import convert_from_path
import optparse

parser = optparse.OptionParser()
parser.add_option('-i', '--input', help='tree file', type='str')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()

drawing = svg2rlg(options.input)
renderPDF.drawToFile(drawing, options.output)
#renderPM.drawToFile(drawing,  options.output + ".png", fmt="PNG")


# images = convert_from_path(options.output)
# for i in range(len(images)):
# 		images[i].save(options.output + str(i) +'.png', 'PNG')
