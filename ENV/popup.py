import sys
import tkinter as tk
import tkinter.font
import tkinter.scrolledtext as tkst

def show_popup(updates):

	window = tk.Tk()
	window.title('Updates summary')
	window.wm_title("Updates summary")
	window.resizable(0,0)

	logpad = tkst.ScrolledText(window,width=50,height=8)
	logpad.grid(row=0,column=0)
	try:	textfile = open(updates,'r')
	except:	
		logpad.delete('1.0',tk.END)
		return
	contents = textfile.read()
	logpad.delete('1.0',tk.END)
	logpad.insert('1.0',contents)
	textfile.close()
	logpad.config(state='disabled')
	window.mainloop()
	return

show_popup(sys.argv[1])
