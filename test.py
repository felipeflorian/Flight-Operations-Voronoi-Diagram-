from HW_2 import *

test_ = FlightOperations("borders_CO.dat", "airports_CO.dat")

# Point 1
test_.plot_vor_airports()

# Point 2
min_, max_ = test_.max_min_areas(True)
print("The airports that reports minimum"
      "and maximum area coverage respectively are: "
      "{0} and {1}".format(min_, max_))

# Point 3
center, r = test_.build_airport(True)
print("The center of the largest circle"
      " is {0} and its radius {1}".format(center, r))

# Point 4
least, most = test_.most_less_crowded(True)
print("The least and most crowded airports"
      " respectively are: {0} and {1}".format(least, most))

# Point 5
a_1, a_2, ch = test_.merge_airport(True)
print("The airports that may be deleted are: {0} and {1}."
      " They will be replaced by an airport with coordinates"
      "{2}".format(a_1, a_2, ch))
