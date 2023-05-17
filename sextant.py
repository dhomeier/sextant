# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Version 0.2 by Derek Homeier <dhomeie@gwdg.de>,
# based on `coordinates-with-pyscript` by Erik Tollerud <erik.tollerud@gmail.com>

# this is needed for conversions that need IERS data
import pyodide_http
import pyodide
pyodide_http.patch_all()

from astropy import coordinates
from astropy.time import Time
from astropy.coordinates.name_resolve import NameResolveError
from astropy.utils.exceptions import AstropyWarning
from astropy.utils.iers.iers import IERSWarning
from erfa.core import ErfaWarning

import warnings

def populate_dropdowns(elem):
    # Hardcode some star names as connecting to CDS seems problematic under PyScript
    star_names = {'Aldebaran': '04h35m55.2390696s +16d30m33.4884636s',
                  'Hamal': '02h07m10.405704s +23d27m44.70318s',
                  'Regulus': '10h08m22.31098512s +11d58m01.9515936s'}

    select_options = [f'<option value="{pos}">{nm}</option>' for nm, pos in star_names.items()]
    
    elem.innerHTML = select_options


datein = Element("date")
timein = Element("time")
lonin = Element("geolon")
latin = Element("geolat")
starin = Element("starpos")
mout = Element("mooncoo")
mdis = Element("moondis")
ebutton = Element("ephembutton")
dbutton = Element("distbutton")
refstar = Element("refstar")
starcoo = Element("starcoo")

populate_dropdowns(refstar.element)

ebutton.element.disabled = False
dbutton.element.disabled = False


def ephem():
    try:
        obsloc = coordinates.EarthLocation.from_geodetic(lonin.element.value, latin.element.value)
    except Exception as e:
        mout.write(f"Fehler bei Positionseingabe: {e} -\t"
                   f"Länge, Breite: {lonin.element.value} {latin.element.value}")
    else:
        mout.write(f"Koordinaten: {obsloc.lon}, {obsloc.lat}")

    with warnings.catch_warnings():
        # Filter warnings about dates out of range (ERFA outside of range 1900-2100)
        warnings.simplefilter("ignore", category=ErfaWarning)
        try:
            obstime = Time(f"{datein.element.value} {timein.element.value}")
        except Exception as e:
            mout.write(f"Fehler bei Datumseingabe: {e} \n{datein.element.value} {timein.element.value}")
        else:
            mout.write(f"Datum, Zeit: {datein.element.value} {timein.element.value} {obstime}")

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=IERSWarning)
                warnings.filterwarnings("ignore", message=".*IERS [dt]a.*", category=AstropyWarning)
                mooncoo = coordinates.get_moon(obstime, obsloc, ephemeris='builtin')
        except Exception as e:
            mout.write(f"Fehler in Mondposition: {e} - Datum, Position: "
                       f"{datein.element.value} {timein.element.value}; {lonin} {latin}")
        else:
            mout.write(f"Mondposition von {obsloc.lon:.3f} Länge {obsloc.lat:+.3f} Breite am "
                   f"{obstime} (UTC): {mooncoo.to_string('hmsdms', precision=2)}")

    return mooncoo

def distanz(mooncoo):
    if len(starin.element.value.strip()) > 2:
        starcoo = None
        try:
            starstr = starin.element.value
            if len(starstr.split()) == 2 and starstr.split()[1].find("d") > 0:
                starcoo = coordinates.SkyCoord(starstr)
            else:
                starcoo = coordinates.SkyCoord.from_name(starstr)
        except NameResolveError:
                mdis.write(f'Stern "{starstr}" ist nicht bekannt.')
        except pyodide.ffi.JsException as e:  # NetworkError?
                mdis.write(f'Keine Verbindung zu CDS für Suche nach "{starstr}" möglich: {e}')
        except Exception as e:
                mdis.write(f"Fehler in Sternposition für {starstr}: {e}")
        else:
            mdis.write(f"Sternposition {starcoo.to_string('hmsdms', precision=2)}")
    else:
        starcoo = coordinates.SkyCoord(refstar.element.value)

    if starcoo is not None:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=ErfaWarning)
                moondis = mooncoo.separation(starcoo)
            mdis.write(f"Abstand zu Sternposition {starcoo.to_string('hmsdms', precision=1)}: "
                       f"    {moondis.to_string(precision=1)}")
        except Exception as e:
            mdis.write(f"Fehler bei Distanzbestimmung zu {starcoo.to_string('hmsdms')}: {e}")
