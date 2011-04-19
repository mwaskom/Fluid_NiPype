import os
import re
import shutil

class HTMLReport(object):

    def __init__(self, path):
        
        self.dir = os.path.dirname(path)
        utildir = os.path.join(self.dir, ".util")
        if not os.path.exists(self.dir):
            os.mkdir(self.dir)
        if not os.path.exists(utildir):
            os.mkdir(utildir)
        self.copy_css()

        self._fid = open(path, "w")

        self._fid.write("<HTML><link REL=stylesheet TYPE=text/css href=.util/report.css>\n")

    def __del__(self):

        self._fid.write("\n</HTML>")
        self._fid.close()

    def copy_css(self):
        
        homedir = os.path.dirname(__file__)
        for util in ["report.css", "seal.png"]:
            shutil.copy(os.path.join(homedir, util),
                        os.path.join(self.dir, ".util", util))
        
    def write_title(self, text):

        self._fid.write("\n<h1>%s</h1>"%text)
                        
    def write_section_head(self, text):

        self._fid.write("\n<h2>%s</h2>"%text)

    def write_text(self, text, bold=False):
       
        bon = ""
        boff = ""
        if bold:
            bon = "<b>"
            boff = "</b>"
        self._fid.write("\n%s%s%s"%(bon,text,boff))

    def newline(self, n=1):

        for l in range(n):
            self._fid.write("\n<br>")

    def write_link(self, text, page, cssclass="body"):

        self._fid.write("\n<p class=%s><a href=%s>%s</a></p>"%(cssclass, page, text))

    def write_link_table(self, link_template, title_template, values, width):

        self._fid.write("\n<p class=link_table>\n")
        for i, val in enumerate(values):
            if not bool(i%width) and i:
                self.newline()
            elif i:
                self._fid.write("&nbsp; -- &nbsp;")
            self._fid.write("\n<a href=%s>%s</a>"%(link_template%val, title_template%val))
        self._fid.write("</p")
    
    def write_image(self, image, size=None):

        if size is None:
            self._fid.write("<br>\n<a href=%s><IMG BORDER=0 SRC=%s ></a><br>\n"%(image, image))
        else:
            self._fid.write("<br>\n"
                            "<a href=%s><IMG BORDER=0 SRC=%s WIDTH=\"%d\" HEIGHT=\"%d\" ></a><br>\n"%(image,image,size[0],size[1]))

    def write_images_across(self, images, size=None):
        
        self.newline()
        for image in images:
            if size is None:
                self._fid.write("\n<a href=%s><IMG BORDER=0 SRC=%s ></a>\n"%(image, image))
            else:
                self._fid.write("\n"
                                "<a href=%s><IMG BORDER=0 SRC=%s WIDTH=\"%d\" HEIGHT=\"%d\" ></a>\n"%(image,image,size[0],size[1]))
        self.newline()
