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


