
# You can also put your plugins in ~/.phy/plugins/.

from phy import IPlugin

# Plugin example:
#
# class MyPlugin(IPlugin):
#     def attach_to_cli(self, cli):
#         # you can create phy subcommands here with click
#         pass

c = get_config()
c.Plugins.dirs = [r'C:\Users\sa4257\.phy\plugins']

# Activate plugins in the Template GUI:
#  - CARFilterPlugin    : adds 'car' / 'car_highpass' raw-data filters (alt+r).
#  - PairCCGPlugin      : registers the PairCorrelogramView (ACG | CCG+jitter |
#                         ACG triptych for exactly 2 selected clusters). Add it
#                         from the GUI View menu as 'PairCorrelogramView'.
c.TemplateGUI.plugins = ['CARFilterPlugin', 'PairCCGPlugin']
