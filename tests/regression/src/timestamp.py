from datetime import datetime

t1 = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
with open("timestamp", "w") as file:
    file.write(str(t1) + "\n")