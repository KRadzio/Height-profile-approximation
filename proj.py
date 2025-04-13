import csv
from matplotlib import pyplot as plt

# ten plik jest do stworzenia wykresów

nodesNum = [5 ,9, 17, 33, 65]
data = {
    "MountEverest.csv": "MountEverest",
    "tczew_starogard.txt": "TczewStarogard",
    "ulm_lugano.txt": "UlmLugano",
    "SpacerniakGdansk.csv": "SpacerniakGdansk"
}
extensions = ["Lagrange", "nodesValues", "Splines"]

# format danych: "x y\n"

for dataName, resultName in data.items():
    for number in nodesNum:
            # dane orginalne
            xOrg = []
            yOrg = []
            with open('{}'.format(dataName), 'r') as file:
                skipLine = False
                delim = ' '
                if dataName == "MountEverest.csv" or dataName == "SpacerniakGdansk.csv":
                    skipLine = True
                    delim = ','
                reader = csv.reader(file, delimiter=delim)
                for line in reader:
                    if(skipLine == True):
                        skipLine = False
                        continue
                    xOrg.append(float(line[0]))
                    yOrg.append(float(line[1]))
            # wyniki interpolacji metodą Lagrange
            xIntL = []
            yIntL = []
            with open('./output/{}{}{}.txt'.format(resultName,number,extensions[0]), 'r') as file:
                reader = csv.reader(file, delimiter=' ')
                for line in reader:
                    xIntL.append(float(line[0]))
                    yIntL.append(float(line[1]))     
            #wyniki interpolacji funkcjami sklejanymi
            xIntS = []
            yIntS = []
            with open('./output/{}{}{}.txt'.format(resultName,number,extensions[2]), 'r') as file:
                reader = csv.reader(file, delimiter=' ')
                for line in reader:
                    xIntS.append(float(line[0]))
                    yIntS.append(float(line[1]))         
            # węzły interpolacji
            xNod = []
            yNod = []
            with open('./output/{}{}{}.txt'.format(resultName,number,extensions[1]), 'r') as file:
                reader = csv.reader(file, delimiter=' ')
                for line in reader:
                    xNod.append(float(line[0]))
                    yNod.append(float(line[1]))
            # wykres L
            fig = plt.figure()      
            if(max(yIntL) > 1.5 * max(yOrg) or min(yIntL) < min(yOrg) / 1.5):
                plt.ylim(min(yOrg) / 1.5, max(yOrg) * 1.5)
            plt.title("Metoda {} dla danych {} ; {} węzłów".format(extensions[0], resultName, number))
            plt.xlabel("Odległość od startu w metrach")
            plt.ylabel("Wysokość w metrach")
            plt.scatter(xOrg, yOrg, color='red', marker='.',label='Wartości rzeczywiste')
            plt.scatter(xIntL, yIntL, color='blue', marker='.',label='Wartości interpolowane')
            plt.scatter(xNod, yNod, marker='x', color='black')
            plt.legend(loc='best', borderaxespad=0.)
            plt.savefig("./charts/{}{}{}.png".format(resultName, number, extensions[0]))
            plt.close(fig)
            #wykres S
            fig = plt.figure()
            plt.title("Metoda {} dla danych {} ; {} węzłów".format(extensions[2], resultName, number))
            plt.xlabel("Odległość od startu w metrach")
            plt.ylabel("Wysokość w metrach")
            plt.scatter(xOrg, yOrg, color='red', marker='.',label='Wartości rzeczywiste')
            plt.scatter(xIntS, yIntS, color='blue', marker='.',label='Wartości interpolowane')
            plt.scatter(xNod, yNod, marker='x', color='black')
            plt.legend(loc='best', borderaxespad=0.)
            plt.savefig("./charts/{}{}{}.png".format(resultName, number, extensions[2]))
            plt.close(fig)