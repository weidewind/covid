from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM
import optparse

parser = optparse.OptionParser()
parser.add_option('-i', '--input', help='tree file', type='str')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()

drawing = svg2rlg(options.input)
renderPDF.drawToFile(drawing, options.output)