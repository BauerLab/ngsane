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
            var n= t.toString().split(".");
            // add comma to decimal positions only
            n[0] = n[0].replace(/\B(?=(\d{3})+(?!\d))/g, ",");
            //Combine both sections again
            return n.join(".");
        }
    } else
        return n;
}

Array.prototype.max = function() {
  return Math.max.apply(null, this);
};

Array.prototype.min = function() {
  return Math.min.apply(null, this);
};

Math.median = function() {
    var ary, numA, i, median;
    ary = Array.prototype.slice.call(arguments);
    for (i = ary.length-1; i >= 0; i--) {if (ary[i] !== +ary[i]) ary[i] = Number.NEGATIVE_INFINITY;}
    numA = function(a, b){return (a-b);};
    ary.sort(numA);
    while (ary.length > 1 && !isFinite(ary[0])) ary.shift();
    median = ary[Math.floor(ary.length/2)];
    if (!isNaN(median))
        return median;
    return "";
}

var isArray = function (obj) {
	return Object.prototype.toString.call(obj) === "[object Array]";
};

Math.mean = function( ){
    var ary, i, mean, sum=0;
    ary = Array.prototype.slice.call(arguments);
    for (i = ary.length-1; i >= 0; i--) {sum += ary[i];}
	mean = (sum / ary.length );
    
    if (!isNaN(mean))
        return mean;
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
    if (!isNaN(v))
        return v;
    return "";
	
}

Math.std = function( ){
    var ary, numA, i;
    ary = Array.prototype.slice.call(arguments);
	var sum = 0;
    for (i = ary.length-1; i >= 0; i--) {sum += ary[i];}
	var avg = (sum / ary.length );
	var i = ary.length;
	var v=0;
	while( i-- ){
		v += Math.pow( (ary[ i ] - avg), 2 );
	}
	return Math.sqrt(v / ary.length);
	
}
