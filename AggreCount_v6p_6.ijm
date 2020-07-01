//------------------------------------------------------------------------
//                    Aggrecount 2020
//------------------------------------------------------------------------
/*
 * Welcome to the Aggrecount Macro developed by the Raman Lab at Tufts University, Boston, MA
 * This macro is intended to quantify immunofluorescent images of cultured cells stained for 
 * the nucleus (eg dapi), aggregates (eg ubiquitin), and cell body (eg Cell Mask). Each image should
 * be a multi-channel or stacked image with each of these stains. 
 * 
 * This macro may be modified to fit your needs by editing the processing steps below. The processing
 * may use any Fiji plugin or function that has been adapted for use in macros. The processing must end
 * with a binary image for nuclei, aggregates and the cell thresholding method. Cell segmentation may 
 * take 8-bit or 16-bit images and requires that the Fiji function "Find Maxima..." is run with the 
 * "segmentation" option be selected.
 */

suffix = "nd2";
 
function AggregateROI(ubq,dapiroi,aggroi,setup,lower,upper,aggcutmin,aggcutmax)
{
	selectWindow(ubq);
	run("Select All");
	run("Duplicate...", "title=template");
	run("16-bit");
	
	//------------------------------------------------------------------------
	//                     AGGREGATE PROCESSING
	//------------------------------------------------------------------------
	//Subtract background as a function of whole image intensity mean

	getStatistics(area, mean, min, max, std, histo);
	run("Subtract...","value="+(area));

	//Enhance contrast
	run("Enhance Contrast...", "saturated="+0.2);

	//rolling ball background subtraction
	run("Subtract Background...", "rolling="+50);

	//Gamma enhancement, values > 1 enhance bright pixels and decrease moderate-low intensity pixels
	//run("Gamma...", "value="+1.34);
	
	//Apply blur
	run("Gaussian Blur...", "sigma=2"); 
	//------------------------------------------------------------------------
	//                     END OF PROCESSING
	//------------------------------------------------------------------------
	
	run("8-bit");

	//Allow manual manipulation of threshold during threshold
	if(setup)
	{
		run("Grays");
		run("Threshold...");
		setThreshold(lower,upper);
		aggnum = 0;
		waitForUser("Click okay when satisfied");
	}
	//Automatic threshold after it has been set
	else 
	{
		setThreshold(lower,upper,"black&white");
	
		run("Convert to Mask");
		run("Analyze Particles...", "size=0-infinity pixel circularity=0.00-1.00 add");
		close("template");

		//Remove aggregates above or below cutoff criteria
		for(j=0;j<roiManager("count");j++)
		{
			roiManager("Select",j);
			getStatistics(area);
			if(area<aggcutmin || area>aggcutmax)
			{
				roiManager("Delete");
				j = j - 1;
			}
		}

		//Get number of aggregate ROI's and save them
		aggnum = roiManager("Count");
		if(aggnum>0) {roiManager("save", aggroi);}
		roiManager("reset");
	}
	
	return aggnum;
}	

function NucleiROI(dapi,dapiroi,nuccut,nucstrictness,batch)
{
	selectWindow(dapi);
	run("Select All");
	run("Duplicate...", "title=template");
	//------------------------------------------------------------------------
	//               NUCLEI PROCESSING
	//------------------------------------------------------------------------
	
	//Enhance contrast
	run("Enhance Contrast...", "saturated=10 normalize");
	
	//Subtracting background as a proportion of the mean
	getStatistics(area, mean, min, max, std, histogram);
	run("Subtract...","value="+(mean/(11-nucstrictness)));
	
	//Median filter to blur and preserve edges
	run("Median...", "radius=13");
	
	//Threshold and turn into binary image
	run("Make Binary");
	run("Dilate");
	run("Fill Holes");
	
	//------------------------------------------------------------------------
	//                  END OF PROCESSING
	//------------------------------------------------------------------------
	
	run("Analyze Particles...", "size=0-infinity pixel circularity=0.00-1.00 add");
	close("template");

	//Remove ROI's that do not meet the cutoff criteria
	selectWindow(dapi);
	for(j=0;j<roiManager("count");j++)
	{
		roiManager("Select",j);
		getStatistics(area);
		if(area<nuccut)
		{
			roiManager("Delete");
			j = j - 1;
		}
	}

	//Display nuclei ROI's if batch mode is not selected
	selectWindow(dapi);
	if(batch==false)
	{
		roiManager("Show All without labels");
		roiManager("Measure");
		waitForUser("Nuclei ROIs");
		run("Clear Results");
	}

	//Get average area of nuclei ROI \
	roiManager("select",Array.getSequence(roiManager("count")));
	if(roiManager("count")>1) {roiManager("COMBINE");}
	getStatistics(area);
	nucarea = area/roiManager("Count");

	//Save ROI's temporarily
	roiManager("save", dapiroi);
	roiManager("reset")
	
	return nucarea;
}

function cellROI(cell,dapiroi,aggroi,cellroi,cellcut,nucarea,batch,segmentation,strictness)
{
	selectWindow(cell);
	run("Select All");
	run("Duplicate...","title=template");

	if(segmentation) //For cell segmentation
	{
		//------------------------------------------------------------------------
		//                     CELL SEGMENTATION
		//------------------------------------------------------------------------
	
		//create a composite image of nuclei and cell bodies
		run("Merge Channels...", "c1="+dapi+" c2="+cell+" create keep");
		run("Stack to RGB");
		run("16-bit");
	
		//Enhance contrast
		run("Enhance Local Contrast (CLAHE)", "blocksize=249 histogram=239 maximum=18 mask=*None* fast_(less_accurate)");
		//Subtract background with rolling ball algorithm
		run("Subtract Background...", "rolling=300");
		//Apply gaussian filter to blur
		run("Gaussian Blur...", "sigma=10"); 
		//Find local maxima and apply pixel intensity watershed algorithm
		run("Find Maxima...", "prominence="+(strictness*4)+" strict exclude output=[Segmented Particles]");
		
		//------------------------------------------------------------------------
		//                     END CELL SEGMENTATION
		//------------------------------------------------------------------------
	}
	else
	{
		//------------------------------------------------------------------------
		//                     CELL THRESHOLDING
		//------------------------------------------------------------------------
		//Enhance contrast
		run("Enhance Contrast...", "saturated=10 normalize");
		//Subtract proportion of mean pixel intensity
		getStatistics(area, mean);
		run("Subtract...","value="+(mean*(strictness/5)));
		//Apply median filter to blur and preserve edges
		run("Median...", "radius=20");
	
		//Turn into binary image and fill in nuclei ROIs
		run("Make Binary");
		roiManager("Open",dapiroi);
		roiManager("Select",Array.getSequence(roiManager("Count")));
		Roi.setFillColor("white");
		roiManager("Fill");
		roiManager("Reset");
	
		//Turn into binary image and dilate
		run("Make Binary");
		run("Dilate");
	
		//------------------------------------------------------------------------
		//                     END CELL THRESHOLDING
		//------------------------------------------------------------------------
	}

	run("Analyze Particles...", "size=0-inifinty pixel circularity=0.00-1.00 add");
	if(segmentation) {close("Composite*");}
	close("template");

	
	//Deletes any cell ROI without a nuclei ROI
	selectWindow(cell);
	rawcellnum = roiManager("count");	
	roiManager("Open",dapiroi);
	totalroinum = roiManager("count");
	roi=0;	
	for(j=0;j<rawcellnum;j++)
	{	
		roiManager("select",roi);
		getStatistics(cellarea, cellmean, min, max, cellstd);
		if(cellarea<cellcut)
		{
			roiManager("Select",roi);
			roiManager("Delete");
			roi=roi-1;
		}
		else 
		{
			containsnuc = false;
			for(k=rawcellnum;k<totalroinum;k++)
			{
				tempdapinum = roiManager("Count") - (totalroinum-k);
				roiManager("select",newArray(tempdapinum,roi));
				roiManager("AND");
				if(selectionType>-1)
				{
					getStatistics(narea, nucmean, min, max, nucstd);
					if(narea>nucarea/5) {containsnuc = true;}
				}
			}
			roiManager("deselect");
			if(containsnuc==false)
			{
				roiManager("Select",roi);
				roiManager("Delete");
				roi = roi - 1;
			}
		}
		roi = roi + 1;
	}

	//Removes nuclei ROI from ROI manager
	roiManager("Select",Array.slice(Array.getSequence(roiManager("Count")),roi,roiManager("Count")));
	roiManager("Delete");

	//Displays cell ROI's if batch mode is not selected
	if(batch==false)
	{
		roiManager("Show All with labels");
		run("Enhance Contrast...", "saturated=0.3");
		roiManager("Measure");
		run("Summarize");
		waitForUser("Cell ROIs");
		run("Clear Results");
	}

	//Saves cell ROI's
	roiManager("save", cellroi);
	roiManager("reset");
	

	return roi;
}

//run("Set Measurements...", "area mean standard max perimeter redirect=None decimal=3");
//Sets options for "Make Binary" function
run("Options...", "iterations=2 count=1 black do=Nothing");

//Define global variables
dir = getDirectory("Choose a folder with images for analysis");
filelist = getFileList(dir);

aggroi = dir+"\\aggroi.zip";
dapiroi = dir+"\\nucroi.zip";
cellroi = dir+"\\cellroi.zip";

rowset = 0;
rowagg = 0;
rowcell = 0;
lower = 60;
upper = 255;
lower1 = 60;
upper1 = 255;
dapislice=newArray(1,1,1);
aggslice=newArray(2,1,1);
cellslice=newArray(3,1,1);

//Begin setup sequence
if(getBoolean("Choose settings for analysis","Setup","Manual input"))
{	
	//Define channel, slice, frame for nuclei, aggregates and cells
	filename=0;
	for(j=0;j<lengthOf(filelist);j++)
	{
		if((endsWith(filelist[j],"tif") || endsWith(filelist[j], suffix)))
		{
			filename = filelist[j];
			j = lengthOf(filelist);			
		}
	}
	if(filename==0) {exit("No images found");}
	
	//Open first image
	bioformsettings = "open=[" + dir + filename + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
	run("Bio-Formats", bioformsettings);

	waitForUser("Select nuclei channel, slice, and frame");
	Stack.getPosition(dapislice[0],dapislice[1],dapislice[2]);
	waitForUser("Select aggregate channel, slice, and frame");
	Stack.getPosition(aggslice[0],aggslice[1],aggslice[2]);
	waitForUser("Select cell body channel, slice, and frame");
	Stack.getPosition(cellslice[0],cellslice[1],cellslice[2]);
	close();

	//Iterate through images to set threshold limits
	while(getBoolean("Set thresholding for aggregates\nCurrent threshold: upper - "+lower+" lower - "+upper+"\nGet image?","Get image","Next step"))
	{
		tempfilepath = File.openDialog("Image for aggregate thresholding");
		bioformsettings = "open=[" + tempfilepath + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
		run("Bio-Formats", bioformsettings);

		Stack.setPosition(aggslice[0],aggslice[1],aggslice[2]);
		aggnum = AggregateROI(getTitle(),dapiroi,aggroi,true,lower,upper,0,1000000);
		getThreshold(lower,upper);
		close("Threshold");
		close("template");
		close("*");
	}
	//Iterate through images to display aggregate sizes and distances from nucleus
	while(getBoolean("Determine aggregate size and distance cut off\nLoad image to measure aggregate size and distance from nucleus","Get image","Next step"))
	{
		tempfilepath = File.openDialog("Image for aggregate size and distance");
		bioformsettings = "open=[" + tempfilepath + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
		run("Bio-Formats", bioformsettings);

		Stack.setPosition(dapislice[0],dapislice[1],dapislice[2]);
		run("Duplicate...","title=template1");
		nucarea = NucleiROI("template1",dapiroi,50,5,true);
		close("template1");
		
		Stack.setPosition(aggslice[0],aggslice[1],aggslice[2]);
		run("Duplicate...","title=template1");
		aggnum = AggregateROI("template1",dapiroi,aggroi,false,lower,upper,0,1000000);
		close("template1");

		roiManager("open",aggroi);
		roiManager("Measure");
		roiManager("open",dapiroi);

		for(j=0;j<aggnum;j++)
		{
			mindis = 10000;
			roiManager("Select",j);
			getSelectionCoordinates(x,y);
			Array.getStatistics(x, min, max, xc, stdDev);
			Array.getStatistics(y, min, max, yc, stdDev);
			
			for(k=aggnum;k<roiManager("Count");k++)
			{
				roiManager("Select",k);
				getSelectionCoordinates(X,Y);
				Array.getStatistics(X, min, max, Xc, stdDev);
				Array.getStatistics(Y, min, max, Yc, stdDev);
		
				if(sqrt(pow(xc-Xc,2) + pow(yc-Yc,2))<500)
				{
					for(l=0;l<X.length;l++)
					{
						for(m=0;m<x.length;m++)
						{
							dis = sqrt(pow(x[m]-X[l],2) + pow(y[m]-Y[l],2));
							if(dis<mindis){mindis=dis;}				
						}
					}
				}
			}
			roiManager("Select",j);
			roiManager("Rename","dis "+mindis+" size "+getResult("Area",j));
		}
	
		run("Make Composite");
		run("From ROI Manager");
		run("Labels...", "color=white font=10 show use draw");
		waitForUser("Check distance and size");
		roiManager("reset");
		run("Clear Results");
		close("*");
	}

	//Present user with strictness options for segementation and thresholding cell bodies
	while(getBoolean("Optimize cell body processing strictness\nThree levels of strictness will be shown: 2.5, 5 & 7.5","Get image","Next step (settings)"))
	{
		tempfilepath = File.openDialog("Image for cell body procesing");
		bioformsettings = "open=[" + tempfilepath + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
		run("Bio-Formats", bioformsettings);
		id=getImageID();
		
		Stack.setPosition(dapislice[0],dapislice[1],dapislice[2]);
		run("Duplicate...","title=Dapi");
		nucarea = NucleiROI("Dapi",dapiroi,50,5,true);
		selectImage(id);
		Stack.setPosition(cellslice[0],cellslice[1],cellslice[2]);
		run("Duplicate...","title=Strictness");
		for(i=2.5;i<10;i=i+2.5)
		{
			//Segmentation
			run("Merge Channels...", "c1=Dapi c2=Strictness create keep");
			run("Stack to RGB");
			run("16-bit");
			run("Enhance Local Contrast (CLAHE)", "blocksize=249 histogram=240 maximum=18 mask=*None* fast_(less_accurate)");
			run("Subtract Background...", "rolling=300");
			run("Median...", "radius=13");
			run("Find Maxima...", "prominence="+(i*4)+" strict exclude output=[Segmented Particles]");
			run("Analyze Particles...", "size=0-inifinty pixel circularity=0.00-1.00 add");
			close("Composite*");
			selectWindow("Strictness");
			run("Duplicate...","title=Strictness_"+i+"_Segmentation");
			run("From ROI Manager");
			run("Overlay Options...", "stroke=none width=0 fill=none");
			roiManager("Reset");

			//Thresholding
			selectWindow("Strictness");
			run("Duplicate...","title=temp");
			run("Enhance Local Contrast (CLAHE)", "blocksize=249 histogram=240 maximum=22 mask=*None* fast_(less_accurate)");
			getStatistics(area, mean);
			run("Subtract...","value="+(mean*(i/5)));
			run("Median...", "radius=13");
			run("Make Binary");
			roiManager("Open",dapiroi);
			roiManager("Select",Array.getSequence(roiManager("Count")));
			Roi.setFillColor("white");
			roiManager("Fill");
			roiManager("Reset");
			run("Make Binary");
			run("Dilate");
			run("Analyze Particles...", "size=0-inifinty pixel circularity=0.00-1.00 add");
			close("temp");
			selectWindow("Strictness");
			run("Duplicate...","title=Strictness_"+i+"_Thresholding");
			if(roiManager("Count")>0)
			{
				run("From ROI Manager");
				run("Overlay Options...", "stroke=none width=0 fill=none");
				roiManager("Reset");
			}

			selectWindow("Strictness");
		}
		selectImage(id);
		close();
		close("Dapi");
		run("Tile");
		waitForUser("Note strictness level from image title");
		close("Strictness*");
	}
	run("Clear Results");
}

if(isOpen("Dataset_summary")) {close("Dataset_summary");}
Table.create("Dataset_summary");
if(isOpen("Dataset_cells")) {close("Dataset_cells");}
Table.create("Dataset_cells");
if(isOpen("Dataset_aggregates")) {close("Dataset_aggregates");}
Table.create("Dataset_aggregates");

//Create the main settings window
Dialog.create("AggreCount settings");
Dialog.setInsets(-10, 120, 3);
Dialog.addMessage("                                       ____                       __   \n"+
"     /\\                               /   ___|                      |   |  \n"+
"    /   \\   _ _  _ _ _ __ _|   |       __  _  _ __   |   |_ \n"+
"   / /\\  \\ /  _  |/ _  |  __/ _ \\  |      / _ \\ | |  |  |  _ \\|  __|\n"+
"  / /_\\  \\  (_| | (_| | | |  __/   |__| (_)|  |_| | | |   |   |_ \n"+
" /_/    \\_\\ _, |\\_, |_|  \\__|\\____\\__/ \\_,_|_| |_|\\__|\n"+
"            _/  | _/  |                                     \n"+
"           |__/ |__/");
Dialog.addMessage("Aggregate thresholding:\n");
Dialog.addNumber("Lower: ",lower);
Dialog.addToSameRow();
Dialog.addNumber("Upper: ",upper);
Dialog.addMessage("Perinuclear distance cutoff\n");
Dialog.addNumber("Distance: ", 10, 0, 4, "pixels");
Dialog.addNumber("Aggresome min size:", 5, 1, 4, "um^2");
Dialog.addNumber("Aggregate min size:", 0.05, 2, 4, "um^2");
Dialog.addToSameRow();
Dialog.addNumber("Max size:", 25, 0, 2, "um^2");
Dialog.addNumber("Nuclei size:", 50, 0, 4, "um^2");
Dialog.addToSameRow();
Dialog.addNumber("Strictness:", 5, 0, 2, "(1-10)");
Dialog.addNumber("Cell size", 75, 0, 4, "um^2");
Dialog.addToSameRow();
Dialog.addNumber("Strictness:", 5, 0, 2, "(1-10)");
Dialog.setInsets(-5,250,3);
Dialog.addMessage("(applies to cell thresholding)");
Dialog.addNumber("Dapi: ",dapislice[0], 0, 2, "channel");
Dialog.addToSameRow();
Dialog.addNumber("",dapislice[1], 0, 2, "slice");
Dialog.addToSameRow();
Dialog.addNumber("",dapislice[2], 0, 2, "frame");
Dialog.addNumber("Aggregates: ",aggslice[0], 0, 2, "channel");
Dialog.addToSameRow();
Dialog.addNumber("",aggslice[1], 0, 2, "slice");
Dialog.addToSameRow();
Dialog.addNumber("",aggslice[2], 0, 2, "frame");
Dialog.addNumber("Cell bodies: ",cellslice[0], 0, 2, "channel");
Dialog.addToSameRow();
Dialog.addNumber("",cellslice[1], 0, 2, "slice");
Dialog.addToSameRow();
Dialog.addNumber("",cellslice[2], 0, 2, "frame");
Dialog.addCheckbox("Find cell bodies?", 1)
Dialog.addChoice("Cell selection method", newArray("Segmentation","Thresholding"));
Dialog.addCheckbox("Save results?",1);
Dialog.addToSameRow();
Dialog.addCheckbox("Batch mode?",0);
Dialog.show();
lower = Dialog.getNumber();
upper = Dialog.getNumber();
perinucdis = Dialog.getNumber();
aggresomecut = Dialog.getNumber();
aggcutmin = Dialog.getNumber();
aggcutmax = Dialog.getNumber();
nuccut = Dialog.getNumber();
nucstrictness = Dialog.getNumber();
cellcut = Dialog.getNumber();
strictness = Dialog.getNumber();
dapislice[0] = Dialog.getNumber();
dapislice[1] = Dialog.getNumber();
dapislice[2] = Dialog.getNumber();
aggslice[0] = Dialog.getNumber();
aggslice[1] = Dialog.getNumber();
aggslice[2] = Dialog.getNumber();
cellslice[0] = Dialog.getNumber();
cellslice[1] = Dialog.getNumber();
cellslice[2] = Dialog.getNumber();
cellfind = Dialog.getCheckbox();
segmentation = Dialog.getChoice();
saveR = Dialog.getCheckbox();
batch = Dialog.getCheckbox();

if(matches(segmentation,"Segmentation")) {segmentation=true;}

if(saveR)
{
	i=0;
	do
	{
		i=i+1;
		newdir = dir+"AC_analysis"+i;
	} while(File.exists(newdir));
	
	File.makeDirectory(newdir);

	setvar = newArray("Lower","Upper","perinuclear distance","aggresome size","minimum aggregate size","maximum aggregate size","minimum nuclei size","nuclei strictness","minimum cell size","cell strictness","nuclei channel","nuclei stack","nuclei frame","aggregate channel","aggregate stack","aggregate frame","cell channel","cell stack","cell frame","Find cells?","Segmentation?","Batch mode");
	setval = newArray(lower,upper,perinucdis,aggresomecut,aggcutmin,aggcutmax,nuccut,nucstrictness,cellcut,strictness,dapislice[0],dapislice[1],dapislice[2],aggslice[0],aggslice[1],aggslice[2],cellslice[0],cellslice[1],cellslice[2],cellfind,segmentation,batch);
	
	Table.create("Settings");
	for(i=0;i<lengthOf(setval);i++)
	{
	Table.set("Variable",i,setvar[i]);
	Table.set("Value",i,setval[i]);
	}
	Table.save(newdir+"\\aggresettings.txt");
	close("Settings");
}

ubq = "Aggregates";
dapi = "Dapi";
cell = "Cells";


//-------------------------------------------------------------------------
//              Image analysis section
//-------------------------------------------------------------------------


//Cycle through each file in directory
for(i=0;i<lengthOf(filelist);i++)
{
filename = filelist[i];
if(saveR) {savepath = newdir+"\\"+substring(filename,0,lengthOf(filename)-4);}

if(batch) {cont = true;}
else      {cont = getBoolean("Next file:\n"+filename+"\nContinue?","Continue","Skip");}  

//Only accept files that end in ND2 or TIF regardless of user input above
if((endsWith(filename,"tif") || endsWith(filename, suffix)) && cont)
{

//Maintain batch mode for ROI manager
close("ROI Manager");
setBatchMode(batch);
roiManager("reset");

//Open next image
bioformsettings = "open=[" + dir + filename + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
run("Bio-Formats", bioformsettings);
id = getImageID();

//opens each channel as a separate image
selectImage(id);
Stack.setPosition(dapislice[0],dapislice[1],dapislice[2]);
run("Duplicate...", "title="+dapi);
selectImage(id);
Stack.setPosition(aggslice[0],aggslice[1],aggslice[2]);
run("Duplicate...", "title="+ubq);
selectImage(id);
if(cellfind)
{
	Stack.setPosition(cellslice[0],cellslice[1],cellslice[2]);
	run("Duplicate...", "title="+cell);
	selectImage(id);
	close();
}
//Determine nuclei ROI (see function above)
nucarea = NucleiROI(dapi,dapiroi,nuccut,nucstrictness,batch);

//Determine aggregate ROI (see function above)
aggnum = AggregateROI(ubq,dapiroi,aggroi,false,lower,upper,aggcutmin,aggcutmax);

//open aggregate and nuclei ROI to determine distances
if(aggnum>0) {roiManager("open",aggroi);}
roiManager("open",dapiroi);
selectWindow(ubq);

aggArea = newArray(aggnum);
aggMean = newArray(aggnum);
aggSTD = newArray(aggnum);
aggLoc = newArray(aggnum);
aggdis = newArray(aggnum);
aggClass = newArray(aggnum);
aggCell = newArray(aggnum);
//Cycle through each aggregate to determine distance from nearest nuclei
for(j=0;j<aggnum;j++)
{
	nuclear = false;
	aggdis[j] = 10000;
	mindisC = 10000;
	centroiddis = 500;
	
	roiManager("Select",j);
	getStatistics(areaj);
	getSelectionCoordinates(x,y);
	Array.getStatistics(x, min, max, xc); //Get x centroid of aggregate
	Array.getStatistics(y, min, max, yc); //get y centroid of aggregate
	
	for(k=aggnum;k<roiManager("Count");k++) //Cycle through each nucleus ROI
	{	
		roiManager("Select",k);
		getSelectionCoordinates(X,Y);
		Array.getStatistics(X, min, max, Xc, stdDev); //get x centroid of nucleus
		Array.getStatistics(Y, min, max, Yc, stdDev); //get y centroid of nucleus
		disC = sqrt(pow(xc-Xc,2) + pow(yc-Yc,2));
		if(disC<mindisC) {mindisC=disC;disCi=k;}
				
		if(disC<centroiddis) //Only calculate distances for agg + nuclei with centroids within 500 pixels
		{
			roiManager("Select",newArray(j,k));
			roiManager("AND");
			getStatistics(area);
			 
			if(area==areaj) //If the aggregate is completely contained within the nuclei
			{
				nuclear = true;
				k = roiManager("Count");
			}
			
			for(l=0;l<X.length;l++) //Nested for loops calc the distance between each point of the agg and nuclei
			{
				for(m=0;m<x.length;m++)
				{
					dis = sqrt(pow(x[m]-X[l],2) + pow(y[m]-Y[l],2));
					if(dis<aggdis[j])
					{aggdis[j]=dis;}
					if(dis==0)
					{k = roiManager("Count");}
				}
			}
		}
		if(k==roiManager("Count") && aggdis[j]==10000) {k=disCi-1;centroiddis=mindisC+1;}
	}
	//adds aggregate data to the arrays
	roiManager("select",j);
	getStatistics(aggArea[j], aggMean[j], min, max, aggSTD[j]);
	if(aggdis[j]<=perinucdis && aggArea[j]>=aggresomecut) {aggClass[j] = "aggresome";}
	else                                            {aggClass[j] = "aggregate";}
	if(nuclear && aggdis[j]>perinucdis)             {aggLoc[j] = "nuclear";}
	else if(aggdis[j]<=perinucdis)                  {aggLoc[j] = "perinuclear";}
	else                                            {aggLoc[j] = "cytosolic";}
}

roiManager("reset");

if(cellfind)
{
	cellnum = cellROI(cell,dapiroi,aggroi,cellroi,cellcut,nucarea,batch,segmentation,strictness);

	roiManager("Open",cellroi);
	
	//Adds cell data to the data table and labels cells
	cellArea = newArray(cellnum);
	cellMean = newArray(cellnum);
	cellSTD = newArray(cellnum);
	for(j=0;j<cellnum;j++)
	{
		roiManager("select",j);
		getStatistics(cellArea[j], cellMean[j], min, max, cellSTD[j]);
	}
	roiManager("reset");

//Beginning of sequence that assigns nuclei and aggregates to cells
roiManager("open", cellroi);
roiManager("open", dapiroi);
nuccellnum = roiManager("Count");
if(aggnum>0) {roiManager("open", aggroi);}

cellnuccount = newArray(cellnum);
cellaggcount = newArray(cellnum);
cellaggarea = newArray(cellnum);
cellcytagg = newArray(cellnum);
cellcytaggarea = newArray(cellnum);
cellpnagg = newArray(cellnum);
cellpnaggarea = newArray(cellnum);
cellnucagg = newArray(cellnum);
cellnucaggarea = newArray(cellnum);
cellaggresomecount = newArray(cellnum);
cellaggresomearea = newArray(cellnum);

//Image variables
cellcount=0;
totagg=0;
totaggarea=0;
totcytagg=0;
totcytaggarea=0;
totpnagg=0;
totpnaggarea=0;
totnucagg=0;
totnucaggarea=0;
totaggresome=0;
totaggresomearea=0;
percentagg=0;
percentaggresome=0;

for(j=0;j<cellnum;j++) //For each cell
{
	for(k=cellnum;k<roiManager("Count");k++) //For each ROI that's NOT a cell
	{
		roiManager('select',newArray(j,k));
		roiManager("AND");
		if(selectionType>-1) //If nuclei or aggregate is contained in the cell
		{
			getStatistics(area);
			if (k<roiManager("count")-aggnum) //If ROI is a nucleus
			{				
				if(area>nucarea){cellnuccount[j] = cellnuccount[j] + floor(area/nucarea);} //If the area is LARGER than an average nucleus, count it as multiple based on size
				else            {cellnuccount[j] = cellnuccount[j] + 1;}
			}
			else  //If ROI is an aggregate
			{
				aggrow = k-nuccellnum;
				
				//Set counters for the cell based on aggregate characteristics
				if (matches(aggClass[aggrow],"aggresome"))
				{
					cellaggresomecount[j] = cellaggresomecount[j] + 1;
					cellaggresomearea[j] = cellaggresomearea[j] + aggArea[aggrow];
				}
				if (matches(aggLoc[aggrow],"perinuclear"))
				{
					cellpnagg[j] = cellpnagg[j] + 1;
					cellpnaggarea[j] = cellpnaggarea[j] + aggArea[aggrow];
				}
				else if (matches(aggLoc[aggrow], "cytosolic"))
				{
					cellcytagg[j] = cellcytagg[j] + 1;
					cellcytaggarea[j] = cellcytaggarea[j] + aggArea[aggrow];
				}
				else 
				{
					cellnucagg[j] = cellnucagg[j] + 1;
					cellnucaggarea[j] = cellnucaggarea[j] + aggArea[aggrow];
				}
				cellaggcount[j] = cellaggcount[j] + 1;
				cellaggarea[j] = cellaggarea[j] + aggArea[aggrow];

				//Assign aggregate to that cell
				aggCell[aggrow] = j;				
			}			
		}
	}

	//Update image variables with new counts
	totagg = totagg + cellaggcount[j];
	totaggarea= totaggarea + cellaggarea[j];
	totaggresome = totaggresome + cellaggresomecount[j];
	totaggresomearea = totaggresomearea + cellaggresomearea[j];
	totcytagg = totcytagg + cellcytagg[j];
	totcytaggarea = totcytaggarea + cellcytaggarea[j];
	totnucagg = totnucagg + cellnucagg[j];
	totnucaggarea = totnucaggarea + cellnucaggarea[j];
	totpnagg = totpnagg + cellpnagg[j];
	totpnaggarea = totpnaggarea + cellpnaggarea[j];
	if(cellaggcount[j]>0 && cellaggcount[j]<cellnuccount[j]){percentagg = percentagg + cellaggcount[j];}
	else if(cellaggcount[j]>0 && cellaggcount[j]>=cellnuccount[j]){percentagg = percentagg + cellnuccount[j];}
	if(cellaggresomecount[j]>0 && cellaggresomecount[j]<cellnuccount[j]){percentaggresome = percentaggresome + cellaggresomecount[j];}
	else if(cellaggresomecount[j]>0 && cellaggresomecount[j]>=cellnuccount[j]){percentaggresome = percentaggresome + cellnuccount[j];}
	cellcount = cellcount + cellnuccount[j];
}


//After cycling through all cells, add image summary data to dataset table

selectWindow("Dataset_summary");
Table.set("File",rowset,filename);
Table.set("Cells",rowset,cellcount);
Table.set("Aggregates per cell",rowset,(totagg/cellcount));
Table.set("Aggregate area per cell",rowset,(totaggarea/cellcount));
Table.set("Cells with aggregates",rowset,((percentagg)));
Table.set("% cells with aggregates",rowset,(percentagg/cellcount));
Table.set("Cells with aggresome",rowset,((percentaggresome)));
Table.set("% cells with aggresome",rowset,(percentaggresome/cellcount));
Table.set("Cells w aggregates but no aggresome",rowset,(percentagg-percentaggresome));
Table.set("% cells w aggregates but no aggresome",rowset,((percentagg-percentaggresome)/cellcount));
Table.set("Average aggregate size",rowset,(totaggarea/totagg));
Table.set("Average aggresome size",rowset,(totaggresomearea/totaggresome));
Table.set("% aggs perinuclear",rowset,(totpnagg/totagg));
Table.set("% aggs cytosolic",rowset,(totcytagg/totagg));
Table.set("% aggs nuclear",rowset,(totnucagg/totagg));
Table.set("Average perinuclear agg size",rowset,(totpnaggarea/totpnagg));
Table.set("Average cytosolic agg size",rowset,(totcytaggarea/totcytagg));
Table.set("Average nuclear agg size",rowset,(totnucaggarea/totnucagg));

Table.update;

//print data to cells window
selectWindow("Dataset_cells");
for(j=0;j<cellnum;j++)
{
	Table.set("File", rowcell, filename);
	Table.set("Nuclei", rowcell, cellnuccount[j]);
	Table.set("Aggregates", rowcell, cellaggcount[j]);
	Table.set("Aggregate area", rowcell, cellaggarea[j]);
	Table.set("Cytosolic aggregates", rowcell, cellcytagg[j]);
	Table.set("Cytosolic agg area", rowcell, cellcytaggarea[j]);
	Table.set("Perinuclear aggregates", rowcell, cellpnagg[j]);
	Table.set("Perinuclear agg area", rowcell, cellpnaggarea[j]);
	Table.set("Nuclear aggregates", rowcell, cellnucagg[j]);
	Table.set("Nuclear agg area", rowcell, cellnucaggarea[j]);
	Table.set("Area", rowcell, cellArea[j]);
	Table.set("Mean", rowcell, cellMean[j]);
	Table.set("StdDev", rowcell, cellSTD[j]);
	rowcell=rowcell+1;
}
Table.update;

rowset = rowset+1; //Proceed to next row in dataset table for next image
}
else //If not aggregate ROIs are present
{
	roiManager("open", dapiroi);
	cellcount = 0;
	for(j=0;j<roiManager("Count");j++)
	{
		roiManager("select",j);
		getStatistics(area);
		tempcount = round(area/nucarea);
		if(tempcount>0) {cellcount = cellcount + tempcount;}
		else {cellcount = cellcount + 1;}
	}
	roiManager("reset");

	totagg=0;
	totaggarea=0;
	totcytagg=0;
	totcytaggarea=0;
	totpnagg=0;
	totpnaggarea=0;
	totnucagg=0;
	totnucaggarea=0;
	totaggresome=0;
	totaggresomearea=0;
	for(j=0;j<lengthOf(aggArea);j++)
	{
		class = aggClass[j];
		loc = aggLoc[j];
		area = aggArea[j];
		if(matches(class,"aggresome"))
		{
			totaggresome = totaggresome + 1;
			totaggresomearea = totaggresomearea + area;
		}
		
		if(matches(loc,"cytosolic"))
		{
			totcytagg = totcytagg + 1;
			totcytaggarea = totcytaggarea + area;
		}
		else if(matches(loc,"perinuclear"))
		{
			totpnagg = totpnagg + 1;
			totpnaggarea = totpnaggarea + area;
		}
		else
		{
			totnucagg = totnucagg + 1;
			totnucaggarea = totnucaggarea + area;
		}

		totagg = totagg + 1;
		totaggarea = totaggarea + area;
	}

	//print data to summary window
	selectWindow("Dataset_summary");
	Table.set("File",rowset,filename);
	Table.set("Cells",rowset,cellcount);
	Table.set("Agg per cell",rowset,(totagg/cellcount));
	Table.set("Agg area per cell",rowset,(totaggarea/cellcount));
	Table.set("Total aggresomes",rowset,(totaggresome));
	Table.set("% cells with aggresome",rowset,(totaggresome/cellcount));
	Table.set("Average aggregate size",rowset,(totaggarea/totagg));
	Table.set("Average aggresome size",rowset,(totaggresomearea/totaggresome));
	Table.set("% aggs perinuclear",rowset,(totpnagg/totagg));
	Table.set("% aggs cytosolic",rowset,(totcytagg/totagg));
	Table.set("% aggs nuclear",rowset,(totnucagg/totagg));
	Table.set("Average perinuclear agg size",rowset,(totpnaggarea/totpnagg));
	Table.set("Average cytosolic agg size",rowset,(totcytaggarea/totcytagg));
	Table.set("Average nuclear agg size",rowset,(totnucaggarea/totnucagg));

	rowset = rowset+1;
	Table.update;
}
//print data to aggregate window
selectWindow("Dataset_aggregates");
for(j=0;j<lengthOf(aggArea);j++)
{
	Table.set("File", rowagg, filename);
	Table.set("Location", rowagg, aggLoc[j]);
	Table.set("Distance",rowagg, aggdis[j]);
	Table.set("Area", rowagg, aggArea[j]);
	Table.set("Mean", rowagg, aggMean[j]);
	Table.set("StdDev",rowagg, aggSTD[j]);
	rowagg=rowagg+1;
}
Table.update;

row=0;
//print data for the specific image file
Table.create("Data");
if(aggnum>0 && cellfind)
{
	for(j=0;j<lengthOf(cellArea);j++)
	{
		Table.set("Class",row,"Cell");
		Table.set("Cell #",row,j+1);
		Table.set("Location",row,"NA");
		Table.set("Area",row,cellArea[j]);
		Table.set("Mean",row,cellMean[j]);
		Table.set("Stdev",row,cellSTD[j]);
		Table.set("Nuclei",row,cellnuccount[j]);
		Table.set("Aggregates",row,cellaggcount[j]);
		Table.set("Aggregate area",row,cellaggarea[j]);
		Table.set("Cytosolic agg",row,cellcytagg[j]);
		Table.set("Cytosolic agg area",row,cellcytaggarea[j]);
		Table.set("Perinuclear agg",row,cellpnagg[j]);
		Table.set("Perinuclear agg area",row,cellpnaggarea[j]);
		Table.set("Nuclear agg",row,cellnucagg[j]);
		Table.set("Nuclear agg area",row,cellnucaggarea[j]);
		Table.set("Aggresomes",row,cellaggresomecount[j]);
		Table.set("Aggresome area",row,cellaggresomearea[j]);
		row=row+1;
	}
}

for(j=0;j<lengthOf(aggArea);j++)
{
	Table.set("Class",row,aggClass[j]);
	Table.set("Cell #",row,aggCell[j]);
	Table.set("Location",row,aggLoc[j]);
	Table.set("Area",row,aggArea[j]);
	Table.set("Mean",row,aggMean[j]);
	Table.set("Stdev",row,aggSTD[j]);
	Table.set("Distance",row,aggdis[j]);
	row=row+1;
}

if(batch==false)
{
	Table.update;
	waitForUser("Click 'OK' for next image");
}

//save function for image data and ROIs
if(saveR)
{
	selectWindow("Data");
	save(savepath+"_analysis.txt");
	close(substring(filename,0,lengthOf(filename)-4)+"_analysis.txt");
	roisave = savepath + "_rois.zip";
	roiManager("Reset");
	if(cellfind){roiManager("Open",cellroi);}
	roiManager("Open",dapiroi);
	if(aggnum>0){roiManager("Open",aggroi);}
	roiManager("save",roisave);
	roiManager("Reset");
}
close("*");

//delete temporary roi files
d=File.delete(aggroi);
d=File.delete(dapiroi);
d=File.delete(cellroi);
}
}
//save dataset windows
if(saveR)
{
	selectWindow("Dataset_summary");
	save(newdir+"\\dataset_summary.txt");

	selectWindow("Dataset_aggregates");
	save(newdir+"\\dataset_aggregates.txt");

	if(cellfind)
	{
		selectWindow("Dataset_cells");
		save(newdir+"\\dataset_cells.txt");
	}
}

