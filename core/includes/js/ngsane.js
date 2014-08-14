function getHash( url ) {
    var hashPos = url.lastIndexOf ( '#' );
    return url.substring( hashPos + 1 );
}

function isNumber(obj) { 
    return !isNaN(parseFloat(obj));
}

// format number with comma
function formatNumber(n) {
    if (isNumber(n))
        return n.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
    else
        return n;
}


