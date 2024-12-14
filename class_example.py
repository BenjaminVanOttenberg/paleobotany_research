class Person:
    def __init__(self, firstname, lastname, age):
        self.firstname = firstname
        self.lastname = lastname
        self.age = age

    def greet(self, message):
        print(f"{self.firstname} {self.lastname} says {message}!")

    def get_fullname(self):
        return f"{self.firstname} {self.lastname}"

    def __str__(self):
        return f"A person named {self.get_fullname()} aged {self.age}"


person = Person("Josiah", "Wang", 20)  # What do you mean I don't look 20? :D
person.greet("I love Machine Learning!")
print(person.age)
print(person.get_fullname())
print(person)
print(type(person))
print(isinstance(person, Person))
print(isinstance(person, object))