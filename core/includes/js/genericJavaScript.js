//This script onclick shows different tabs
$( "div.panel div.tabs a" ).click(function() {

    //Finds old selected tablink (heading) and sets to unselected
//    console.log($(this).parent().find("a.selected"));
    $(this).parent().find("a.selected").attr('class', "");
    
    //Sets new tab heading to selected
//    console.log($(this).attr('href'));
    $(this).attr('class', "selected");
//    console.log($(this));
    
    //Using id from tab heading find content div id
    var id = getHash ($(this).attr('href'));
    var ContentDivId = "DC_"+id;
    //console.log(ContentDivId);
    
    //Sets old content div to tabContent hide
//    console.log($('#' + ContentDivId).parent().find("div.tabContent"));
    $('#' + ContentDivId).parent().find("div.tabContent").attr('class', "tabContent_hide");
    
    //Sets new content div to tabContent (able to be seen)
    //console.log($('#' + ContentDivId).attr('class'));
    $('#' + ContentDivId).attr('class', "tabContent");
    //console.log($('#' + ContentDivId).attr('class'));

});

// add filtering for invalid rows
$.fn.dataTableExt.afnFiltering.push(
	 function( oSettings, aData, iDataIndex ) {
		
		var search = $("#search");
		if ($("#search")[0].value=="") {
             console.log("No search value");
			return true;
		} else {
			var searchArray = $("#search")[0].value.split(/\s+/);
			console.log(searchArray);
			var sArrayLength = searchArray.length;
			var numColumns = oSettings.aoColumns.length;
			//check relevant columns using 
			// cycle through columns
			var display = 0;
			//document.write(searchArray[0]);
			
			//document.write(sArrayLength);
			//document.write(numColumns);
			var iColumn = 0;
			while (iColumn < numColumns && display == 0) {
						
				//go through search array and check if cell regex matches any search items..
				for (isTerm = 0; isTerm < sArrayLength; isTerm++) {
					//alert (
					//document.write(searchArray[isTerm].length);
					if (searchArray[isTerm].length != 0) {
						var patt1=new RegExp(searchArray[isTerm]);
						//document.write(patt1);
						if(patt1.test(aData[iColumn])) {
							display = 1;
						}
					}

				}
				iColumn++;
			}
			if (display == 1) {
				return true;
			} else {
				return false;
			}
		}
	}
);