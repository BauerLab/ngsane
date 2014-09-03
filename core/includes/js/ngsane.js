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
    var ary, numA, i;
    ary = Array.prototype.slice.call(arguments);
    for (i = ary.length-1; i >= 0; i--) {if (ary[i] !== +ary[i]) ary[i] = Number.NEGATIVE_INFINITY;}
    numA = function(a, b){return (a-b);};
    ary.sort(numA);
    while (ary.length > 1 && !isFinite(ary[0])) ary.shift();
    return ary[Math.floor(ary.length/2)];
}