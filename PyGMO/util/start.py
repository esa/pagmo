class analysis:
    def __init__(self):
        pass    
    def start(self):
        import Tkinter as tk 
        from Tkinter import Tk, Frame, Checkbutton, Label, Entry, Radiobutton
        from Tkinter import IntVar,BooleanVar, StringVar, BOTH


        class parameters(Frame):
  
            def __init__(self, parent):
                Frame.__init__(self, parent)   
                 
                self.parent = parent        
                self.initUI()
                
            def initUI(self):
              
                self.parent.title("ANALYSIS")

                self.pack(fill=BOTH, expand=1)

                self.b1=BooleanVar()#sample
                self.s1=StringVar()#number of points
                self.i1=IntVar()#sobol/lhs

                self.b2=BooleanVar()#plot f-distributions
                self.s2=StringVar()#percentiles

                self.b3=BooleanVar()#perform lin_conv test
                self.s3=StringVar()#number of pairs

                self.b41=BooleanVar()#linear reg
                self.b42=BooleanVar()#linear reg w/ interaction
                self.s42=StringVar()#order of interaction
                self.b43=BooleanVar()#poly reg
                self.s43=StringVar()#degree of regression

                self.b51=BooleanVar()#SVM
                self.s51=StringVar()#threshold svm
                self.b52=BooleanVar()#DAC
                self.s52=StringVar()#threshold dac
                self.s53=StringVar()#ktest svn
                self.s54=StringVar()#ktune svn
                self.s55=StringVar()#ktest dac

                self.b6=BooleanVar()#perform local search
                self.s6=StringVar()#number of initial points
                self.b60=BooleanVar()#cluster results
                self.i61=IntVar()#var_ratio/k
                self.s61=StringVar()#var_ratio
                self.s62=StringVar()#k
                self.b63=BooleanVar()#scatter plot
                self.s63=StringVar()#scatter plot dimensions
                self.b64=BooleanVar()#PCP together
                self.b65=BooleanVar()#PCP separate

                self.b71=BooleanVar()
                self.b72=BooleanVar()


                #SAMPLING

                self.l1=Label(self, text="SAMPLING")
                self.l1.grid(row=0,sticky='W')

                self.cb1 = Checkbutton(self, text="Sample",
                    variable=self.b1,command=self.add1)
                self.cb1.select()
                self.cb1.grid(row=1,column=0,sticky='W')

                self.l11=Label(self, text="N. points:")
                self.l11.grid(row=1,column=0, sticky='E')

                self.e1=Entry(self,textvariable=self.s1,width=5, justify=tk.CENTER)
                self.e1.insert(0,'1000')
                self.e1.grid(row=1,column=1,sticky='W')


                self.r11=Radiobutton(self, text='sobol',variable=self.i1, value=0)
                self.r11.grid(row=1,column=1,sticky='e')

                self.r12=Radiobutton(self, text='lhs',variable=self.i1, value=1)
                self.r12.grid(row=1,column=2,sticky='w')

                #F-DISTRIBUTION

                self.l2=Label(self,text='F-DISTRIBUTION FEATURES')
                self.l2.grid(row=2,sticky='W')

                self.cb2 = Checkbutton(self, text="Plot f-distributions",
                    variable=self.b2)
                self.cb2.grid(row=3,column=0,sticky='W')

                self.l2=Label(self, text="Show percentiles:")
                self.l2.grid(row=3,column=1,sticky='E')

                self.e2=Entry(self,textvariable=self.s2,width=15, justify=tk.CENTER)
                self.e2.insert(0,'0,5,10,25,50,100')
                self.e2.grid(row=3,column=2,sticky='w')

                #DEGREE OF LINEARITY AND CONVEXITY

                self.l3=Label(self,text='LINEARITY AND CONVEXITY')
                self.l3.grid(row=4,sticky='W')

                self.cb3 = Checkbutton(self, text="Perform test",
                    variable=self.b3,command=self.add3)
                self.cb3.grid(row=5,column=0,sticky='W')

                self.l3=Label(self, text="Number of pairs:",state=tk.DISABLED)
                self.l3.grid(row=5,column=1,sticky='E')

                self.e3=Entry(self,textvariable=self.s3,width=5,state=tk.DISABLED, justify=tk.CENTER)
                self.e3.grid(row=5,column=2,sticky='w')

                #META-MODEL FEATURES

                self.l4=Label(self,text='META-MODEL FEATURES')
                self.l4.grid(row=6,sticky='W')

                self.cb41 = Checkbutton(self, text="Linear regression",
                    variable=self.b41)
                self.cb41.grid(row=7,column=0,sticky='W')

                self.cb42 = Checkbutton(self, text="Linear w/ interaction",
                    variable=self.b42,command=self.add42)
                self.cb42.grid(row=8,column=0,sticky='W')

                self.l42=Label(self, text="Order(s) of interaction:",state=tk.DISABLED)
                self.l42.grid(row=8,column=1,sticky='E')

                self.e42=Entry(self,textvariable=self.s42,width=5,state=tk.DISABLED, justify=tk.CENTER)
                self.e42.insert(0,'2')
                self.e42.grid(row=8,column=2,sticky='w')

                self.cb43 = Checkbutton(self, text="Polynomial regression",
                    variable=self.b43,command=self.add43)
                self.cb43.grid(row=9,column=0,sticky='W')

                self.l43=Label(self, text="Degree(s) of regression:",state=tk.DISABLED)
                self.l43.grid(row=9,column=1,sticky='E')

                self.e43=Entry(self,textvariable=self.s43,width=5,state=tk.DISABLED, justify=tk.CENTER)
                self.e43.insert(0,'2')
                self.e43.grid(row=9,column=2,sticky='w')
 

                #MULTI-MODALITY

                self.l5=Label(self, text="LEVELSET FEATURES / MMI")
                self.l5.grid(row=10,sticky='W')

                self.cb51 = Checkbutton(self, text="SVM",
                    variable=self.b51,command=self.add51)
                self.cb51.grid(row=11,column=0,sticky='W')

                self.l51=Label(self, text="Threshold(s):",state=tk.DISABLED)
                self.l51.grid(row=11,column=0,sticky='E')

                self.e51=Entry(self,textvariable=self.s51,state=tk.DISABLED,width=10, justify=tk.CENTER)
                self.e51.grid(row=11,column=1,sticky='w')

                self.cb52 = Checkbutton(self, text="DAC",
                    variable=self.b52,command=self.add52)
                self.cb52.grid(row=12,column=0,sticky='W')

                self.l52=Label(self, text="Threshold(s):",state=tk.DISABLED)
                self.l52.grid(row=12,column=0,sticky='E')

                self.e52=Entry(self,textvariable=self.s52,state=tk.DISABLED,width=10, justify=tk.CENTER)
                self.e52.grid(row=12,column=1,sticky='w')

                self.l53=Label(self, text="K_test:",state=tk.DISABLED)
                self.l53.grid(row=11,column=1,sticky='E')

                self.e53=Entry(self,textvariable=self.s53,state=tk.DISABLED,width=5, justify=tk.CENTER)
                self.e53.grid(row=11,column=2,sticky='w')

                self.l54=Label(self, text="K_tune:",state=tk.DISABLED)
                self.l54.grid(row=11,column=2,sticky='E')

                self.e54=Entry(self,textvariable=self.s54,state=tk.DISABLED,width=5, justify=tk.CENTER)
                self.e54.grid(row=11,column=3,sticky='w')

                self.l55=Label(self, text="K_test:",state=tk.DISABLED)
                self.l55.grid(row=12,column=1,sticky='E')

                self.e55=Entry(self,textvariable=self.s55,state=tk.DISABLED,width=5, justify=tk.CENTER)
                self.e55.grid(row=12,column=2,sticky='w')

                #LOCAL SEARCH
                
                self.l6=Label(self, text="LOCAL SEARCH")
                self.l6.grid(row=13,sticky='W')

                self.cb6 = Checkbutton(self, text="Perform local search",
                    variable=self.b6,command=self.add6)
                self.cb6.grid(row=14,column=0,sticky='W')

                self.l60=Label(self, text="N. initial points:",state=tk.DISABLED)
                self.l60.grid(row=14, column=1,sticky='e')

                self.e6=Entry(self,textvariable=self.s6,state=tk.DISABLED,width=5, justify=tk.CENTER)
                self.e6.grid(row=14,column=2,sticky='w')

                self.cb60 = Checkbutton(self, text="Cluster results",
                    variable=self.b60,command=self.add61,state=tk.DISABLED)
                self.cb60.grid(row=15,column=0,sticky='W')

                self.r61=Radiobutton(self, text="Fix variance ratio:",state=tk.DISABLED,variable=self.i61,value=0,command=self.add62)
                self.r61.grid(row=15,column=1,sticky='w')

                self.e61=Entry(self,textvariable=self.s61,state=tk.DISABLED,width=5, justify=tk.CENTER)
                self.e61.grid(row=15,column=2,sticky='w')

                self.r62=Radiobutton(self, text="Fix number of clusters:",state=tk.DISABLED,variable=self.i61,value=1,command=self.add62)
                self.r62.grid(row=16,column=1,sticky='w')

                self.e62=Entry(self,textvariable=self.s62,state=tk.DISABLED,width=5, justify=tk.CENTER)
                self.e62.grid(row=16,column=2,sticky='w')

                self.cb63 = Checkbutton(self, text="Scatter plot",
                    variable=self.b63,command=self.add63,state=tk.DISABLED)
                self.cb63.grid(row=16,column=3,sticky='W')

                self.l63=Label(self, text="Dim(s):",state=tk.DISABLED)
                self.l63.grid(row=16,column=3,sticky='E')

                self.e63=Entry(self,textvariable=self.s54,state=tk.DISABLED,width=10, justify=tk.CENTER)
                self.e63.grid(row=16,column=4,sticky='w')

                self.cb64 = Checkbutton(self, text="Plot PCP (together)",
                    variable=self.b64,state=tk.DISABLED)
                self.cb64.grid(row=14,column=3,sticky='W')

                self.cb65 = Checkbutton(self, text="Plot PCP (per cluster)",
                    variable=self.b65,state=tk.DISABLED)
                self.cb65.grid(row=15,column=3,sticky='W')

                #CURVATURE

                self.l7=Label(self, text="CURVATURE")
                self.l7.grid(row=17,sticky='W')

                self.l71=Label(self, text="  -Gradient")
                self.l71.grid(row=18, sticky='W')

                self.cb71 = Checkbutton(self, text="Get gradient",
                    variable=self.b71)
                self.cb71.select()
                self.cb71.grid(row=19,column=0,sticky='W')

                self.l72=Label(self, text="  -Hessian")
                self.l72.grid(row=20, sticky='W')

                self.cb72 = Checkbutton(self, text="Get hessian",
                    variable=self.b72)
                self.cb72.select()
                self.cb72.grid(row=21,column=0,sticky='W')

            def add1(self):
                if self.b1.get()==True:
                    self.l11.configure(state=tk.NORMAL)
                    self.e1.configure(state=tk.NORMAL)
                    self.r11.configure(state=tk.NORMAL)
                    self.r12.configure(state=tk.NORMAL)

                else:
                    self.l11.configure(state=tk.DISABLED)
                    self.e1.configure(state=tk.DISABLED)
                    self.r11.configure(state=tk.DISABLED)
                    self.r12.configure(state=tk.DISABLED)

            def add3(self):
                if self.b3.get()==True:
                    self.l3.configure(state=tk.NORMAL)
                    self.e3.configure(state=tk.NORMAL)
                    self.e3.insert(0,'X')

                else:
                    self.l3.configure(state=tk.DISABLED)
                    self.e3.delete(0,500)
                    self.e3.configure(state=tk.DISABLED)

            def add42(self):
                if self.b42.get()==True:
                    self.l42.configure(state=tk.NORMAL)
                    self.e42.configure(state=tk.NORMAL)
                    self.e42.insert(0,'2')

                else:
                    self.l42.configure(state=tk.DISABLED)
                    self.e42.delete(0,500)
                    self.e42.configure(state=tk.DISABLED)


            def add43(self):
                if self.b43.get()==True:
                    self.l43.configure(state=tk.NORMAL)
                    self.e43.configure(state=tk.NORMAL)
                    self.e43.insert(0,'2')

                else:
                    self.l43.configure(state=tk.DISABLED)
                    self.e43.delete(0,500)
                    self.e43.configure(state=tk.DISABLED)

            def add51(self):
                if self.b51.get()==True:
                    self.l51.configure(state=tk.NORMAL)
                    self.e51.configure(state=tk.NORMAL)
                    self.e51.insert(0,'50')
                    self.l53.configure(state=tk.NORMAL)
                    self.e53.configure(state=tk.NORMAL)
                    self.e53.insert(0,'10')
                    self.l54.configure(state=tk.NORMAL)
                    self.e54.configure(state=tk.NORMAL)
                    self.e54.insert(0,'3')
                else:
                    self.l51.configure(state=tk.DISABLED)
                    self.e51.delete(0,500)
                    self.e51.configure(state=tk.DISABLED)
                    self.l53.configure(state=tk.DISABLED)
                    self.e53.delete(0,500)
                    self.e53.configure(state=tk.DISABLED)
                    self.l54.configure(state=tk.DISABLED)
                    self.e54.delete(0,500)                    
                    self.e54.configure(state=tk.DISABLED)


            def add52(self):
                if self.b52.get()==True:
                    self.l52.configure(state=tk.NORMAL)
                    self.e52.configure(state=tk.NORMAL)
                    self.e52.insert(0,'50')
                    self.l55.configure(state=tk.NORMAL)
                    self.e55.configure(state=tk.NORMAL)
                    self.e55.insert(0,'10')
                else:
                    self.l52.configure(state=tk.DISABLED)
                    self.e52.delete(0,500)
                    self.e52.configure(state=tk.DISABLED)
                    self.l55.configure(state=tk.DISABLED)
                    self.e55.delete(0,500)
                    self.e55.configure(state=tk.DISABLED)

            def add6(self):
                if self.b6.get()==True:
                    self.r61.configure(state=tk.NORMAL)
                    self.r61.select()
                    self.r62.configure(state=tk.NORMAL)
                    self.e61.configure(state=tk.NORMAL)
                    self.e61.insert(0,'0.95')
                    self.cb60.configure(state=tk.NORMAL)
                    self.cb60.select()
                    self.e62.delete(0,500)
                    self.e62.configure(state=tk.DISABLED)
                    self.cb63.configure(state=tk.NORMAL)
                    self.cb63.deselect()
                    self.l63.configure(state=tk.DISABLED)
                    self.e63.delete(0,500)
                    self.e63.configure(state=tk.DISABLED)
                    self.cb64.configure(state=tk.NORMAL)
                    self.cb65.configure(state=tk.NORMAL)
                    self.cb64.deselect()
                    self.cb65.deselect()
                    self.l60.configure(state=tk.NORMAL)
                    self.e6.configure(state=tk.NORMAL)
                    self.e6.insert(0,'X')
                else:
                    self.r61.configure(state=tk.DISABLED)
                    self.r62.configure(state=tk.DISABLED)
                    self.e61.delete(0,500)
                    self.e61.configure(state=tk.DISABLED)
                    self.cb60.deselect()
                    self.cb60.configure(state=tk.DISABLED)
                    self.e62.delete(0,500)
                    self.e62.configure(state=tk.DISABLED)
                    self.cb63.deselect()
                    self.cb63.configure(state=tk.DISABLED)
                    self.l63.configure(state=tk.DISABLED)
                    self.e63.delete(0,500)
                    self.e63.configure(state=tk.DISABLED)
                    self.cb64.deselect()
                    self.cb65.deselect()
                    self.cb64.configure(state=tk.DISABLED)
                    self.cb65.configure(state=tk.DISABLED)
                    self.l60.configure(state=tk.DISABLED)
                    self.e6.delete(0,500)
                    self.e6.configure(state=tk.DISABLED)

            def add61(self):
                if self.b60.get()==True:                
                    self.r61.configure(state=tk.NORMAL)
                    self.r61.select()
                    self.r62.configure(state=tk.NORMAL)
                    self.e61.configure(state=tk.NORMAL)
                    self.e61.insert(0,'0.95')
                    self.e62.configure(state=tk.DISABLED)
                    self.cb63.configure(state=tk.NORMAL)
                    self.cb63.deselect()
                    self.l63.configure(state=tk.DISABLED)
                    self.e63.delete(0,500)
                    self.e63.configure(state=tk.DISABLED)
                    self.cb64.configure(state=tk.NORMAL)
                    self.cb65.configure(state=tk.NORMAL)
                    self.cb64.deselect()
                    self.cb65.deselect()
                else:
                    self.r61.configure(state=tk.DISABLED)
                    self.e61.delete(0,500)
                    self.e61.configure(state=tk.DISABLED)
                    self.r62.configure(state=tk.DISABLED)
                    self.e62.delete(0,500)
                    self.e62.configure(state=tk.DISABLED)
                    self.cb63.configure(state=tk.DISABLED)
                    self.cb63.deselect()
                    self.l63.configure(state=tk.DISABLED)
                    self.e63.delete(0,500)
                    self.e63.configure(state=tk.DISABLED)
                    self.cb64.deselect()
                    self.cb65.deselect()
                    self.cb64.configure(state=tk.DISABLED)
                    self.cb65.configure(state=tk.DISABLED)

            def add62(self):
                if self.i61.get()==0:
                    self.e61.configure(state=tk.NORMAL)
                    if self.e61.get()=='':
                        self.e61.insert(0,'0.95')
                    self.e62.delete(0,500)
                    self.e62.configure(state=tk.DISABLED)
                else:
                    self.e62.configure(state=tk.NORMAL)
                    if self.e62.get()=='':
                        self.e62.insert(0,'10')
                    self.e61.delete(0,500)
                    self.e61.configure(state=tk.DISABLED)

            def add63(self):
                if self.b63.get()==True:
                    self.l63.configure(state=tk.NORMAL)
                    self.e63.delete(0,500)
                    self.e63.configure(state=tk.NORMAL)
                else:
                    self.l63.configure(state=tk.DISABLED)
                    self.e63.delete(0,500)
                    self.e63.configure(state=tk.DISABLED)

        root = Tk()
        root.geometry("800x600+300+300")
        app = parameters(root)
        root.mainloop()