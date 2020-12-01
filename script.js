var pi = 3.14159265358979;

var sm_a = 6378137.0;
var sm_b = 6356752.314;
var sm_EccSquared = 6.69437999013e-3;

var UTMScaleFactor = 0.9996;

function DegToRad(deg) {
  return (deg / 180.0) * pi;
}

function RadToDeg(rad) {
  return (rad / pi) * 180.0;
}

function ArcLengthOfMeridian(phi) {
  var alpha, beta, gamma, delta, epsilon, n;
  var result;

  n = (sm_a - sm_b) / (sm_a + sm_b);

  alpha =
    ((sm_a + sm_b) / 2.0) *
    (1.0 + Math.pow(n, 2.0) / 4.0 + Math.pow(n, 4.0) / 64.0);

  beta =
    (-3.0 * n) / 2.0 +
    (9.0 * Math.pow(n, 3.0)) / 16.0 +
    (-3.0 * Math.pow(n, 5.0)) / 32.0;

  gamma = (15.0 * Math.pow(n, 2.0)) / 16.0 + (-15.0 * Math.pow(n, 4.0)) / 32.0;

  delta =
    (-35.0 * Math.pow(n, 3.0)) / 48.0 + (105.0 * Math.pow(n, 5.0)) / 256.0;

  epsilon = (315.0 * Math.pow(n, 4.0)) / 512.0;

  result =
    alpha *
    (phi +
      beta * Math.sin(2.0 * phi) +
      gamma * Math.sin(4.0 * phi) +
      delta * Math.sin(6.0 * phi) +
      epsilon * Math.sin(8.0 * phi));

  return result;
}


function UTMCentralMeridian(zone) {
  var cmeridian;
  cmeridian = DegToRad(-183.0 + zone * 6.0);

  return cmeridian;
}

function FootpointLatitude(y) {
  var y_, alpha_, beta_, gamma_, delta_, epsilon_, n;
  var result;

  n = (sm_a - sm_b) / (sm_a + sm_b);

  alpha_ =
    ((sm_a + sm_b) / 2.0) * (1 + Math.pow(n, 2.0) / 4 + Math.pow(n, 4.0) / 64);

  y_ = y / alpha_;

  beta_ =
    (3.0 * n) / 2.0 +
    (-27.0 * Math.pow(n, 3.0)) / 32.0 +
    (269.0 * Math.pow(n, 5.0)) / 512.0;

  gamma_ = (21.0 * Math.pow(n, 2.0)) / 16.0 + (-55.0 * Math.pow(n, 4.0)) / 32.0;

  delta_ =
    (151.0 * Math.pow(n, 3.0)) / 96.0 + (-417.0 * Math.pow(n, 5.0)) / 128.0;

  epsilon_ = (1097.0 * Math.pow(n, 4.0)) / 512.0;

  result =
    y_ +
    beta_ * Math.sin(2.0 * y_) +
    gamma_ * Math.sin(4.0 * y_) +
    delta_ * Math.sin(6.0 * y_) +
    epsilon_ * Math.sin(8.0 * y_);

  return result;
}


function MapLatLonToXY(phi, lambda, lambda0, xy) {
  var N, nu2, ep2, t, t2, l;
  var l3coef, l4coef, l5coef, l6coef, l7coef, l8coef;
  var tmp;

  ep2 = (Math.pow(sm_a, 2.0) - Math.pow(sm_b, 2.0)) / Math.pow(sm_b, 2.0);

  nu2 = ep2 * Math.pow(Math.cos(phi), 2.0);

  N = Math.pow(sm_a, 2.0) / (sm_b * Math.sqrt(1 + nu2));
  
  t = Math.tan(phi);
  t2 = t * t;
  tmp = t2 * t2 * t2 - Math.pow(t, 6.0);

  l = lambda - lambda0;

  l3coef = 1.0 - t2 + nu2;
  l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
  l5coef = 5.0 - 18.0 * t2 + t2 * t2 + 14.0 * nu2 - 58.0 * t2 * nu2;
  l6coef = 61.0 - 58.0 * t2 + t2 * t2 + 270.0 * nu2 - 330.0 * t2 * nu2;
  l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - t2 * t2 * t2;
  l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - t2 * t2 * t2;

  xy[0] =
    N * Math.cos(phi) * l +
    (N / 6.0) * Math.pow(Math.cos(phi), 3.0) * l3coef * Math.pow(l, 3.0) +
    (N / 120.0) * Math.pow(Math.cos(phi), 5.0) * l5coef * Math.pow(l, 5.0) +
    (N / 5040.0) * Math.pow(Math.cos(phi), 7.0) * l7coef * Math.pow(l, 7.0);

  xy[1] =
    ArcLengthOfMeridian(phi) +
    (t / 2.0) * N * Math.pow(Math.cos(phi), 2.0) * Math.pow(l, 2.0) +
    (t / 24.0) * N * Math.pow(Math.cos(phi), 4.0) * l4coef * Math.pow(l, 4.0) +
    (t / 720.0) * N * Math.pow(Math.cos(phi), 6.0) * l6coef * Math.pow(l, 6.0) +
    (t / 40320.0) *
      N *
      Math.pow(Math.cos(phi), 8.0) *
      l8coef *
      Math.pow(l, 8.0);

  return;
}


function MapXYToLatLon(x, y, lambda0, philambda) {
  var phif, Nf, Nfpow, nuf2, ep2, tf, tf2, tf4, cf;
  var x1frac, x2frac, x3frac, x4frac, x5frac, x6frac, x7frac, x8frac;
  var x2poly, x3poly, x4poly, x5poly, x6poly, x7poly, x8poly;

  phif = FootpointLatitude(y);

  ep2 = (Math.pow(sm_a, 2.0) - Math.pow(sm_b, 2.0)) / Math.pow(sm_b, 2.0);

  cf = Math.cos(phif);

  nuf2 = ep2 * Math.pow(cf, 2.0);

  Nf = Math.pow(sm_a, 2.0) / (sm_b * Math.sqrt(1 + nuf2));
  Nfpow = Nf;

  tf = Math.tan(phif);
  tf2 = tf * tf;
  tf4 = tf2 * tf2;
  
  x1frac = 1.0 / (Nfpow * cf);

  Nfpow *= Nf; 
  x2frac = tf / (2.0 * Nfpow);

  Nfpow *= Nf; 
  x3frac = 1.0 / (6.0 * Nfpow * cf);

  Nfpow *= Nf; 
  x4frac = tf / (24.0 * Nfpow);

  Nfpow *= Nf; 
  x5frac = 1.0 / (120.0 * Nfpow * cf);

  Nfpow *= Nf; 
  x6frac = tf / (720.0 * Nfpow);

  Nfpow *= Nf; 
  x7frac = 1.0 / (5040.0 * Nfpow * cf);

  Nfpow *= Nf; 
  x8frac = tf / (40320.0 * Nfpow);
 
  x2poly = -1.0 - nuf2;

  x3poly = -1.0 - 2 * tf2 - nuf2;

  x4poly =
    5.0 +
    3.0 * tf2 +
    6.0 * nuf2 -
    6.0 * tf2 * nuf2 -
    3.0 * (nuf2 * nuf2) -
    9.0 * tf2 * (nuf2 * nuf2);

  x5poly = 5.0 + 28.0 * tf2 + 24.0 * tf4 + 6.0 * nuf2 + 8.0 * tf2 * nuf2;

  x6poly = -61.0 - 90.0 * tf2 - 45.0 * tf4 - 107.0 * nuf2 + 162.0 * tf2 * nuf2;

  x7poly = -61.0 - 662.0 * tf2 - 1320.0 * tf4 - 720.0 * (tf4 * tf2);

  x8poly = 1385.0 + 3633.0 * tf2 + 4095.0 * tf4 + 1575 * (tf4 * tf2);
  
  philambda[0] =
    phif +
    x2frac * x2poly * (x * x) +
    x4frac * x4poly * Math.pow(x, 4.0) +
    x6frac * x6poly * Math.pow(x, 6.0) +
    x8frac * x8poly * Math.pow(x, 8.0);

  philambda[1] =
    lambda0 +
    x1frac * x +
    x3frac * x3poly * Math.pow(x, 3.0) +
    x5frac * x5poly * Math.pow(x, 5.0) +
    x7frac * x7poly * Math.pow(x, 7.0);

  return;
}


function LatLonToUTMXY(lat, lon, zone, xy) {
  MapLatLonToXY(lat, lon, UTMCentralMeridian(zone), xy);

  xy[0] = xy[0] * UTMScaleFactor + 500000.0;
  xy[1] = xy[1] * UTMScaleFactor;
  if (xy[1] < 0.0) xy[1] = xy[1] + 10000000.0;

  return zone;
}


function UTMXYToLatLon(x, y, zone, southhemi, latlon) {
  var cmeridian;
  x -= 500000.0;
  x /= UTMScaleFactor;
  if (southhemi) y -= 10000000.0;
  y /= UTMScaleFactor;
  cmeridian = UTMCentralMeridian(zone);
  MapXYToLatLon(x, y, cmeridian, latlon);

  return;
}


function btnToUTM_OnClick() {
  var xy = new Array(2);

  if (isNaN(parseFloat(document.getElementById("txtLongitude").value))) {
    alert("Please enter a valid longitude in the lon field.");
    return false;
  }

  lon = parseFloat(document.getElementById("txtLongitude").value);

  if (lon < -180.0 || 180.0 <= lon) {
    alert(
      "The longitude you entered is out of range.  " +
        "Please enter a number in the range [-180, 180)."
    );
    return false;
  }

  if (isNaN(parseFloat(document.getElementById("txtLatitude").value))) {
    alert("Please enter a valid latitude in the lat field.");
    return false;
  }

  lat = parseFloat(document.getElementById("txtLatitude").value);

  if (lat < -90.0 || 90.0 < lat) {
    alert(
      "The latitude you entered is out of range.  " +
        "Please enter a number in the range [-90, 90]."
    );
    return false;
  }

  // Compute the UTM zone.
  zone = Math.floor((lon + 180.0) / 6) + 1;

  zone = LatLonToUTMXY(DegToRad(lat), DegToRad(lon), zone, xy);

  /* Set the output controls.  */
  document.getElementById("txtX").value = xy[0];
  document.getElementById("txtY").value = xy[1];
  document.getElementById("txtZone").value = zone;
  if (lat < 0)
    // Set the S button.
    document.getElementById("rbtnHemisphereS").checked = true;
  // Set the N button.
  else document.getElementById("rbtnHemisphereN").checked = true;

  return true;
}


function btnToGeographic_OnClick() {
  latlon = new Array(2);
  var x, y, zone, southhemi;

  if (isNaN(parseFloat(document.getElementById("txtX").value))) {
    alert("Please enter a valid easting in the x field.");
    return false;
  }

  x = parseFloat(document.getElementById("txtX").value);

  if (isNaN(parseFloat(document.getElementById("txtY").value))) {
    alert("Please enter a valid northing in the y field.");
    return false;
  }

  y = parseFloat(document.getElementById("txtY").value);

  if (isNaN(parseInt(document.getElementById("txtZone").value))) {
    alert("Please enter a valid UTM zone in the zone field.");
    return false;
  }

  zone = parseFloat(document.getElementById("txtZone").value);

  if (zone < 1 || 60 < zone) {
    alert(
      "The UTM zone you entered is out of range.  " +
        "Please enter a number in the range [1, 60]."
    );
    return false;
  }

  if (document.getElementById("rbtnHemisphereS").checked == true) southhemi = true;
  else southhemi = false;

  UTMXYToLatLon(x, y, zone, southhemi, latlon);

  document.getElementById("txtLongitude").value = RadToDeg(latlon[1]);
  document.getElementById("txtLatitude").value = RadToDeg(latlon[0]);

  return true;
}
