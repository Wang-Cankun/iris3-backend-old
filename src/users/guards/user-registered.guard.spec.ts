import { UserRegisteredGuard } from './user-registered.guard';

describe('UserRegisteredGuard', () => {
  it('should be defined', () => {
    expect(new UserRegisteredGuard()).toBeDefined();
  });
});
