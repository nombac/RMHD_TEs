"""
Target configuration for thermal equilibrium plotting.
Converted from targ_hirose2015.pro
"""

# Data path configuration
DATA_BASE_PATH = './data'

# Color definitions
COLORS = {
    'h2006': 'limegreen',
    'h2007': 'red',
    'h2009': 'orange',
    'h2011': 'violet',
    'h2014': 'blue',
    'h2015': 'black',
    'h2016': 'brown',
    'dummy': 'black'
}

# Target definitions
# Each target is a dictionary with:
# - name: simulation name
# - ave: [start, end] time averaging window in orbital periods
# - color_key: key for color lookup
TARGETS = [
    {'name': 'h2006', 'ave': [50, 150], 'color_key': 'h2006'},
    {'name': 'h2007', 'ave': [50, 150], 'color_key': 'h2007'},
    
    # h2009 series
    {'name': '0211b', 'ave': [50, 150], 'color_key': 'h2009'},
    {'name': '0519b', 'ave': [50, 200], 'color_key': 'h2009'},
    {'name': '1126b', 'ave': [100, 400], 'color_key': 'h2009'},
    {'name': '0520a', 'ave': [150, 300], 'color_key': 'h2009'},
    {'name': '0320a', 'ave': [150, 300], 'color_key': 'h2009'},
    {'name': '090304a', 'ave': [100, 400], 'color_key': 'h2009'},
    {'name': '090423a', 'ave': [50, 250], 'color_key': 'h2009'},
    
    # h2016 series (self-gravity)
    {'name': 'self6004', 'ave': [20, 100], 'color_key': 'h2016'},  # Sigma = 80
    {'name': 'self6025', 'ave': [20, 100], 'color_key': 'h2016'},  # Sigma = 90
    {'name': 'self6000', 'ave': [20, 120], 'color_key': 'h2016'},  # Sigma = 100
    {'name': 'self6005', 'ave': [10, 100], 'color_key': 'h2016'},  # Sigma = 150
    {'name': 'self6006', 'ave': [10, 100], 'color_key': 'h2016'},  # Sigma = 200
    {'name': 'self6026', 'ave': [10, 100], 'color_key': 'h2016'},  # Sigma = 230
    {'name': 'self6007', 'ave': [10, 100], 'color_key': 'h2016'},  # Sigma = 250
    
    # h2011 series
    {'name': '200811z', 'ave': [100, 200], 'color_key': 'h2011'},
    
    # h2015 series
    {'name': 'ws0800', 'ave': [100, 200], 'color_key': 'h2015'},
    {'name': 'ws0804', 'ave': [50, 200], 'color_key': 'h2015'},
    {'name': 'ws0806', 'ave': [100, 200], 'color_key': 'h2015'},
    {'name': 'ws0809', 'ave': [20, 120], 'color_key': 'h2015'},
    {'name': 'ws0818', 'ave': [50, 150], 'color_key': 'h2015'},
    {'name': 'ws0819', 'ave': [150, 250], 'color_key': 'h2015'},
    {'name': 'ws0820', 'ave': [10, 110], 'color_key': 'h2015'},
    {'name': 'ws0829', 'ave': [20, 140], 'color_key': 'h2015'},
    {'name': 'ws0834', 'ave': [10, 110], 'color_key': 'h2015'},
    {'name': 'ws0837', 'ave': [10, 125], 'color_key': 'h2015'},
    {'name': 'ws0847', 'ave': [20, 100], 'color_key': 'h2015'},
    {'name': 'ws0848', 'ave': [20, 120], 'color_key': 'h2015'},
    {'name': 'ws0852', 'ave': [10, 110], 'color_key': 'h2015'},
    {'name': 'ws0855', 'ave': [50, 150], 'color_key': 'h2015'},
    {'name': 'ws0856', 'ave': [81, 181], 'color_key': 'h2015'},
    {'name': 'ws0814', 'ave': [10, 100], 'color_key': 'h2015'},
    {'name': 'ws0822', 'ave': [30, 120], 'color_key': 'h2015'},
    {'name': 'ws0805', 'ave': [80, 180], 'color_key': 'h2015'},
    {'name': 'ws0803', 'ave': [60, 160], 'color_key': 'h2015'},
    {'name': 'ws0844', 'ave': [25, 115], 'color_key': 'h2015'},
    {'name': 'ws0812', 'ave': [30, 150], 'color_key': 'h2015'},
    {'name': 'ws0826', 'ave': [40, 140], 'color_key': 'h2015'},
    {'name': 'ws0843', 'ave': [15, 125], 'color_key': 'h2015'},
    {'name': 'wb0828', 'ave': [50, 150], 'color_key': 'h2015'},
    {'name': 'wc0828', 'ave': [50, 200], 'color_key': 'h2015'},
    {'name': 'ws0860', 'ave': [20, 120], 'color_key': 'h2015'},
    {'name': 'ws0802', 'ave': [35, 130], 'color_key': 'h2015'},
    {'name': 'ws0832', 'ave': [15, 115], 'color_key': 'h2015'},
    {'name': 'ws0831', 'ave': [20, 120], 'color_key': 'h2015'},
    {'name': 'ws0859', 'ave': [50, 150], 'color_key': 'h2015'},
    
    # h2014 series (upper branch)
    {'name': 'ws0430', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0429', 'ave': [90, 190], 'color_key': 'h2014'},
    {'name': 'ws0439', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0491', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0470', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0492', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0425', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0427', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0437', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0433', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0436', 'ave': [30, 130], 'color_key': 'h2014'},
    {'name': 'wt0487', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0446', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'wt0442', 'ave': [10, 110], 'color_key': 'h2014'},
    
    # h2014 series (middle branch)
    {'name': 'ws0474', 'ave': [80, 180], 'color_key': 'h2014'},
    {'name': 'ws0473', 'ave': [50, 150], 'color_key': 'h2014'},
    
    # h2014 series (lower branch)
    {'name': 'ws0466', 'ave': [20, 100], 'color_key': 'h2014'},
    {'name': 'ws0438', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'wt0435', 'ave': [120, 220], 'color_key': 'h2014'},
    {'name': 'ws0462', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0434', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0476', 'ave': [50, 150], 'color_key': 'h2014'},
    {'name': 'ws0445', 'ave': [50, 150], 'color_key': 'h2014'},
    
    # Dummy entry to mark end
    {'name': 'dummy', 'ave': [0, 0], 'color_key': 'dummy'}
]
