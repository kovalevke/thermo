import cantera as ct
import math

#функция по расчёту равновесной скорости звука
def equilSoundSpeeds(gas, rtol=1.0e-6, max_iter=5000):
    
    gas.equilibrate('TP', rtol=rtol, max_iter=max_iter)

    s0 = gas.s
    p0 = gas.P
    r0 = gas.density
    
    p1 = p0*1.0001

    gas.SP = s0, p1

    gas.equilibrate('SP', rtol=rtol, max_iter=max_iter)
    aequil = math.sqrt((p1 - p0)/(gas.density - r0))
    
    return aequil

def findPressureInCritical(gas):
    # запоминаем параметры в КС
    s0 = gas.s
    h0 = gas.h
    p0 = gas.P
    
    mdot = 1        # задаём расход газа 1кг/с
    # задаём первое приближение для площади и давления в критике
    amin = 1.e14    # принимаем площадь F'' побольше
    p_crit = p0     # принимаем давление в критике такое же как в камере
    
    for r in range(500,1000):    # делаем разбиение давления на 1000 частей
                                # и считаем давления от 0,5*p_k и выше
        p = p0*(r+1)/1001.0
        # рассчитываем параметры газа при заданном давлении
        gas.SP = s0, p
        gas.equilibrate('SP')
        W2 = 2.0*(h0 - gas.h)      # h + W^2/2 = h0 - закон сохранения энергии, учебник Дорофеева стр.274 (3-изд)
        W = math.sqrt(W2)
        area = mdot/(gas.density*W) # расчёт F'' из ур-ия неразрывности
        
        if (area <= amin):
            amin = area
            p_crit = p
            W_crit = W
        else:
            break
    print(f'Скорость ПС в критике: {W_crit:.2f}')
    gas.SP = s0, p_crit
    gas.equilibrate('SP')
    print(f'Скорость звука в критике: {equilSoundSpeeds(gas):.2f}')
    
    return p_crit

print('====================Расчет топливной пары Метан-Кислород====================')


gas = ct.Solution('nasa_gas_2.yaml')
fuel = "CH4"
oxidizer = "O2"
alpha=0.8
gas.set_equivalence_ratio(1/(alpha), fuel, oxidizer)
k0 = gas.stoich_air_fuel_ratio(fuel, oxidizer, basis='mass')

p_k=1*1e6
#p_a=100000

H_gor =-5566
H_ok =-398.3

km=k0*alpha
print(f'Коэффициент избытка окислителя ={km}')

#Расчёт массовых долей
#    Массовая доля CH4:
mol_gor=(1/(1+km))
print(f'Массовая доля CH4:{mol_gor}')
#    Массовая доля O2:
mol_ok=1*km/(1+km)
print(f'Массовая доля O2:{mol_ok}')

#Общая энтальпия смеси:
H_sum=((mol_gor*H_gor)+(mol_ok*H_ok))*1000
print(f'Общая энтальпия смеси: {H_sum} Дж/кг')



gas.TP = 300, p_k
gas.equilibrate('TP')
gas.HP = H_sum, p_k
gas.equilibrate('HP')
k = gas.cp/gas.cv
print(gas())

#Запоминаем параметры в камере для удобства использования в последующих расчётах
T_k = gas.T
k_k = gas.cp/gas.cv
R_k = gas.cp-gas.cv
V_k = gas.v
S0 = gas.s

print("")
print("=========================Расчёт до критики:=========================")
print('----------Сечение в камере---------- ')

#рассчитываем равновесную скорость звука
A_zv = equilSoundSpeeds(gas)
#так как в процессе расчёта "испортили" исходные параметры, то возвращаем всё как было
gas.SP = S0, p_k
gas.equilibrate('SP')

#Выводим на экран что насчитали
print(f'Давление = {p_k/1000000} МПа')
print(f'Температура = {T_k:.2f} К')
print(f'Показатель адиабаты (зам): {k_k:.3f}')
print(f'Газовая постоянная: {R_k:.2f}')
print(f'Cp: {gas.cp:.2f}')
print(f'Cv: {gas.cv:.2f}')
print(f'Скорость звука: {A_zv:.2f}')
#print(gas.viscosity)
#print(gas.thermal_conductivity)

print("")

#расчёт равновесной Cv
gas.TD = T_k*1.01, 1/V_k
gas.equilibrate('TV')
U2 = gas.int_energy_mass
gas.TD = T_k*0.99, 1/V_k
gas.equilibrate('TV')
U1 = gas.int_energy_mass;
CVEQ = (U2-U1)/(.02*T_k);

gas.SP = S0, p_k
gas.equilibrate('SP')

#расчёт равновесной Cp
gas.TP = T_k*1.01, p_k
gas.equilibrate('TP')
H2 = gas.enthalpy_mass
gas.TP = T_k*0.99, p_k
gas.equilibrate('TP')
H1 = gas.enthalpy_mass
CPEQ = (H2-H1)/(.02*T_k)

#расчёт равновесного показателя адиабаты и вывод всех параметров
kEQ = CPEQ/CVEQ
print(f'Cp_eq: {CPEQ:.2f}')
print(f'Cv_eq: {CVEQ:.2f}')
print(f'Показатель адиабаты (равн): {kEQ:.3f}')

#так как в процессе расчёта "испортили" исходные параметры, то возвращаем всё как было
gas.SP = S0, p_k
gas.equilibrate('SP')

print("")
#рассчитываем давление в критике
p_kp = findPressureInCritical(gas)
print(f'Давление в критике: {p_kp/1000000:.5f}')
#так как в процессе расчёта "испортили" исходные параметры, то возвращаем всё как было
gas.SP = S0, p_k
gas.equilibrate('SP')
