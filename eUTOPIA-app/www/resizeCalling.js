function resizePlot() {
	var plotImg = this.getElementsByTagName('IMG')[0];
	nWidth=this.clientWidth-5;
	nHeight=this.clientHeight-5;
	nWidthStr=nWidth + 'px';
	nHeightStr=nHeight + 'px';
	plotImg.style.width=nWidthStr;
	plotImg.style.height=nHeightStr;
	alert("here")
}
//var plotDiv = document.getElementsByClassName('shiny-plot-output')[0];
//addResizeListener(plotDiv, resizePlot);
//var plotDivs = document.getElementsByClassName('shiny-plot-output');
var plotDivs = document.getElementsByClassName('sizeable');
alert(plotDivs.length)
var i;
for (i=0; i<plotDivs.length; i++) {
    addResizeListener(plotDivs[i], resizePlot);
    //plotDivs[i].onresize = function() {resizePlot()};
}
