from constructor import *
from element import *
from functors import *
from graded_ring_element import *
from hecke_triangle_groups import *
from graded_ring import *
from space import *
from abstract_ring import *
from abstract_space import *

from analytic_type import *
from series_constructor import *

from sage.categories.pushout import pushout, construction_tower
import pdb
cm = sage.structure.element.get_coercion_model()

AT = AnalyticType()
JFC = JFSeriesConstructor()

K = JFC.K()
wp = JFC.wp()

JF = ModularForms()
(x,y,z,d,a,b) = JF._pol_ring.gens()
JFE = JF.extend_type("weak", ring=True)
jinv = JFE(x**JFE._group.n()/(x**JFE._group.n()-y**2))
jinvred = jinv.reduce()
