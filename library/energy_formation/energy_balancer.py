import json

class Energyknockdown(object):
    def energyknockdownretriever(self):
        energy_knockdown_hash = {
        "hb":
          {
          "CAT_1":0.9,
          "CAT_2":0.9,
          "CAT_3":0.9,
          "CAT_4":0.9,
          "CAT_5":0.9,
          "CAT_6":0.9,
          "CAT_7":0.9,
          "CAT_8":0.9,
          "CAT_9":0.9
        },
          "vdw":
          {
          "CAT_1":0.1,
          "CAT_2":0.1,
          "CAT_3":0.1,
          "CAT_4":0.1,
          "CAT_5":0.1,
          "CAT_6":0.1,
          "CAT_7":0.1,
          "CAT_8":0.1,
          "CAT_9":0.1
          }
        }

        print energy_knockdown_hash
        return energy_knockdown_hash
