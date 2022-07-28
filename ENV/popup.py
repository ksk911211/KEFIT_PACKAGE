import sys, os
import tkinter as tk
import tkinter.font
import tkinter.scrolledtext as tkst

def show_popup(updates):

	if not os.path.isfile(updates): return
	textfile = open(updates,'r')
	contents = textfile.read().split('\n')
	textfile.close()
	name = contents[0]
	contents.remove(name)
	content = '';
	for item in contents: content += item+'\n'
	window = tk.Tk()
	window.title('%s Updates'%name)
	window.wm_title("%s Updates"%name)
	window.resizable(0,0)

	logpad = tkst.ScrolledText(window,width=50,height=20)
	logpad.grid(row=0,column=0)
	logpad.delete('1.0',tk.END)
	logpad.insert('1.0',content)
	logpad.config(state='disabled')

	b1 = tk.Button(window, text="OK", bg = "lightgray",command=lambda: close_popup(window),height = 1,width = 5)
	b1.grid(row=1, column=0)
	window.mainloop()

def close_popup(window):
	window.destroy()

show_popup(sys.argv[1])
