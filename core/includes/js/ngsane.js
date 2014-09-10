function getHash( url ) {
    var hashPos = url.lastIndexOf ( '#' );
    return url.substring( hashPos + 1 );
}

// format number with comma
function formatNumber(n) {
    // only attempt conversion if n is numberic
    if (!isNaN(parseFloat(n))){
        var t = Math.abs(parseFloat(n));
        if (t == 0)
            return n;
        else if (t < 0.01)
            return t.toExponential();
        else {
            // treat pre and post differently
            var n= (Math.round(t * 100) / 100).toString().split(".");
            // add comma to decimal positions only
            n[0] = n[0].replace(/\B(?=(\d{3})+(?!\d))/g, ",");
            //Combine both sections again
            return n.join(".");
        }
    } else
        return n;
}

Math.max = function() {
    var ary, m;
    ary = Array.prototype.slice.call(arguments);
    numA = function(a, b){return (b-a);};
    ary.sort(numA);
    m = ary[0];
    if (!isNaN(m) && isFinite(m))
        return m;
    return "";
};

Math.min = function() {
    var ary, m;
    ary = Array.prototype.slice.call(arguments);
    numA = function(a, b){return (a-b);};
    ary.sort(numA);
    m = ary[0];
    if (!isNaN(m) && isFinite(m))
        return m;
    return "";
};

Math.median = function() {
    var ary, numA, i, m;
    ary = Array.prototype.slice.call(arguments);
    numA = function(a, b){return (a-b);};
    ary.sort(numA);
    m = ary[Math.floor(ary.length/2)];
    if (!isNaN(m) && isFinite(m))
        return m;
    return "";
}


var isArray = function (obj) {
	return Object.prototype.toString.call(obj) === "[object Array]";
};

Math.mean = function( ){
    var ary, i, m, sum=0;
    ary = Array.prototype.slice.call(arguments);
    for (i = ary.length-1; i >= 0; i--) {sum += ary[i];}
	m = (sum / ary.length );
    
    if (!isNaN(m) && isFinite(m))
        return m;
    return "";

}

Math.variance = function( ){
    var ary, i, mean, v=0, sum=0;
    ary = Array.prototype.slice.call(arguments);
    for (i = ary.length-1; i >= 0; i--) {sum += ary[i];}
	mean = (sum / ary.length );
    console.log(ary)
    console.log(mean)
    for (i = ary.length-1; i >= 0; i--) {
		v += Math.pow( (ary[ i ] - mean), 2 );
	}
	v /= ary.length;
    if (!isNaN(v) && isFinite(v))
        return v;
    return "";
	
}

Math.std = function( ){
    var ary, avg, i, sddev, v=0, sum=0 ;
    ary = Array.prototype.slice.call(arguments);
    for (i = ary.length-1; i >= 0; i--) {sum += ary[i];}
	avg = (sum / ary.length );
    i = ary.length;
	while( i-- ){
		v += Math.pow( (ary[ i ] - avg), 2 );
	}
	sddev = Math.sqrt(v / ary.length);
    if (!isNaN(sddev) && isFinite(sddev))
        return sddev;
    return "";
	
}
